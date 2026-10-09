(* eigencube.ml - Minimalistic Rubik's Cube Solver in OCaml
   A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
   and multi-phase A* search with move-budgeted restarts. *)

open Base
open Stdio
open Poly

type vec = int * int * int
type mat = vec * vec * vec
type move = { normal : vec; dir : int }

type cube =
  mat array (* 26 cubelets, each mapped to its current 3x3 rotation matrix *)

let norm1 (x, y, z) = Int.abs x + Int.abs y + Int.abs z
let dot (x1, y1, z1) (x2, y2, z2) = (x1 * x2) + (y1 * y2) + (z1 * z2)

let ( *@ ) ((m00, m01, m02), (m10, m11, m12), (m20, m21, m22)) (x, y, z) =
  ( (m00 * x) + (m01 * y) + (m02 * z),
    (m10 * x) + (m11 * y) + (m12 * z),
    (m20 * x) + (m21 * y) + (m22 * z)
  )

let transpose ((a, b, c), (d, e, f), (g, h, i)) =
  ((a, d, g), (b, e, h), (c, f, i))

let ( *@* ) (r0, r1, r2) m =
  let c0, c1, c2 = transpose m in
  ( (dot r0 c0, dot r0 c1, dot r0 c2),
    (dot r1 c0, dot r1 c1, dot r1 c2),
    (dot r2 c0, dot r2 c1, dot r2 c2)
  )

let diag (x, y, z) = ((x, 0, 0), (0, y, 0), (0, 0, z))
let id3 : mat = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
let crange = [ -1; 0; 1 ]

let all_vectors =
  List.concat_map crange ~f:(fun x ->
      List.concat_map crange ~f:(fun y ->
          List.map crange ~f:(fun z -> (x, y, z))
      )
  )

let cubelets =
  all_vectors |> List.filter ~f:(fun v -> norm1 v > 0) |> Array.of_list

let num_cubelets = Array.length cubelets
let unit_vectors = List.filter all_vectors ~f:(fun v -> norm1 v = 1)

let moves =
  List.concat_map unit_vectors ~f:(fun normal ->
      [ -1; 1 ] |> List.map ~f:(fun dir -> { normal; dir })
  )
  |> Array.of_list

let num_moves = Array.length moves

let rot_mat { normal = x, y, _; dir } =
  if x <> 0 then ((1, 0, 0), (0, 0, dir), (0, -dir, 0))
  else if y <> 0 then ((0, 0, dir), (0, 1, 0), (-dir, 0, 0))
  else ((0, dir, 0), (-dir, 0, 0), (0, 0, 1))

let rot_matrices = Array.map moves ~f:rot_mat

let inv_move =
  Array.init num_moves ~f:(fun m ->
      let inv = { normal = moves.(m).normal; dir = -moves.(m).dir } in
      let rec find i = if moves.(i) = inv then i else find (i + 1) in
      find 0
  )

let opposite_pruned lm m =
  let x1, y1, z1 = moves.(lm).normal and x2, y2, z2 = moves.(m).normal in
  (x1 > x2 || (x1 = x2 && (y1 > y2 || (y1 = y2 && z1 > z2))))
  && x1 = -x2
  && y1 = -y2
  && z1 = -z2

let is_cubelet_solved c r =
  let colors = diag c in
  let sticker_directions = r *@* colors in
  sticker_directions = colors

let is_cubelet_pos_solved c r = r *@ c = c
let solved_cube () : cube = Array.create ~len:num_cubelets id3

let apply_move m (cube : cube) : cube =
  let v = moves.(m).normal and rm = rot_matrices.(m) in
  Array.init num_cubelets ~f:(fun i ->
      let r = cube.(i) in
      if dot v (r *@ cubelets.(i)) > 0 then rm *@* r else r
  )

let is_cube_solved (cube : cube) =
  let rec loop i =
    i = num_cubelets || (is_cubelet_solved cubelets.(i) cube.(i) && loop (i + 1))
  in
  loop 0

let hash_mat ((a, b, c), (d, e, f), (g, h, i)) =
  a
  + (3 * b)
  + (9 * c)
  + (27 * d)
  + (81 * e)
  + (243 * f)
  + (729 * g)
  + (2187 * h)
  + (6561 * i)

let hash_cube (arr : cube) =
  let h = ref 17 in
  for i = 0 to num_cubelets - 1 do
    h := (!h * 31) + hash_mat arr.(i)
  done;
  !h

module CubeTbl = Stdlib.Hashtbl.Make (struct
  type t = cube

  let equal = Poly.equal
  let hash = hash_cube
end)

(* Lazy distance heuristics *)
let dist_solved_cache = Hashtbl.Poly.create ()
let dist_pos_cache = Hashtbl.Poly.create ()

let single_cubelet_bfs c r is_goal =
  if is_goal c r then 0
  else
    let q = Queue.create () in
    let visited = Hashtbl.Poly.create () in
    Queue.enqueue q (r, 0);
    Hashtbl.set visited ~key:r ~data:0;
    let res = ref None in
    while Option.is_none !res && not (Queue.is_empty q) do
      let curr_r, d = Queue.dequeue_exn q in
      for m = 0 to num_moves - 1 do
        if Option.is_none !res then begin
          let next_r = rot_matrices.(m) *@* curr_r in
          if not (Hashtbl.mem visited next_r) then
            begin if is_goal c next_r then res := Some (d + 1)
            else begin
              Hashtbl.set visited ~key:next_r ~data:(d + 1);
              Queue.enqueue q (next_r, d + 1)
            end
            end
        end
      done
    done;
    match !res with
    | Some d -> d
    | None -> 0

let min_moves_to_solved c r =
  Hashtbl.find_or_add dist_solved_cache (c, r) ~default:(fun () ->
      single_cubelet_bfs c r is_cubelet_solved
  )

let min_moves_to_pos c r =
  Hashtbl.find_or_add dist_pos_cache (c, r) ~default:(fun () ->
      single_cubelet_bfs c r is_cubelet_pos_solved
  )

let top_layer_heuristic (cube : cube) =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let _, _, z = cubelets.(i) in
    if z = 1 then
      let d = min_moves_to_solved cubelets.(i) cube.(i) in
      sum := !sum +. (Float.of_int d **. p)
  done;
  (!sum **. (1.0 /. p)) /. 8.0

let middle_layer_heuristic (cube : cube) =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let _, _, z = cubelets.(i) in
    if z >= 0 then
      let d = min_moves_to_solved cubelets.(i) cube.(i) in
      sum := !sum +. (Float.of_int d **. p)
  done;
  (!sum **. (1.0 /. p)) /. 4.0

let bottom_layer_edge_heuristic (cube : cube) =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    let _, _, z = c in
    if not (z = -1 && norm1 c = 3) then
      let d = min_moves_to_solved c cube.(i) in
      sum := !sum +. (Float.of_int d **. p)
  done;
  (!sum **. (1.0 /. p)) /. 3.0

let bottom_layer_corner_heuristic (cube : cube) =
  let p = 0.5 and s1 = ref 0.0 and s2 = ref 0.0 and s3 = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) and r = cube.(i) in
    let _, _, z = c in
    if z = 1 then s1 := !s1 +. (Float.of_int (min_moves_to_solved c r) **. p)
    else if z = 0 then
      s2 := !s2 +. (Float.of_int (min_moves_to_solved c r) **. p)
    else if z = -1 then
      let d =
        if norm1 c = 3 then min_moves_to_pos c r else min_moves_to_solved c r
      in
      s3 := !s3 +. (Float.of_int d **. p)
  done;
  ((!s1 **. (1.0 /. p)) /. 5.0)
  +. ((!s2 **. (1.0 /. p)) /. 3.0)
  +. ((!s3 **. (1.0 /. p)) /. 8.0)

type 'a heap = Empty | Node of float * 'a * 'a heap list

let empty_heap = Empty

let merge_heap h1 h2 =
  match (h1, h2) with
  | Empty, h | h, Empty -> h
  | Node (p1, x1, l1), Node (p2, x2, l2) ->
    if Float.(p1 <= p2) then Node (p1, x1, h2 :: l1) else Node (p2, x2, h1 :: l2)

let push_heap h p x = merge_heap (Node (p, x, [])) h

let rec merge_pairs = function
  | [] -> Empty
  | [ h ] -> h
  | h1 :: h2 :: hs -> merge_heap (merge_heap h1 h2) (merge_pairs hs)

let pop_heap = function
  | Empty -> None
  | Node (_, x, hs) -> Some (x, merge_pairs hs)

let random_gauss mean std =
  let u1 = Float.max 1e-15 (Stdlib.Random.float 1.0)
  and u2 = Stdlib.Random.float 1.0 in
  mean
  +. std
     *. (Float.sqrt (-2.0 *. Float.log u1) *. Float.cos (2.0 *. Float.pi *. u2))

let total_moves_simulated = ref 0

let reconstruct came_from dst =
  let rec loop curr acc =
    match CubeTbl.find_opt came_from curr with
    | Some (p, mv) -> loop p (mv :: acc)
    | None -> acc
  in
  (dst, loop dst [])

let should_prune last_move m =
  match last_move with
  | None -> false
  | Some lm -> m = inv_move.(lm) || opposite_pruned lm m

let astar start is_goal heuristic random_weight max_moves =
  if is_goal start then Some (start, [])
  else
    let rec attempt budget =
      let frontier = ref (push_heap empty_heap 0.0 start) in
      let came_from = CubeTbl.create 8192 in
      let cost_so_far = CubeTbl.create 8192 in
      CubeTbl.replace cost_so_far start 0;
      let simulated = ref 0 in
      let frontier_exhausted = ref false in
      let solution = ref None in

      let step_move src last_move m =
        if not (should_prune last_move m) then begin
          Int.incr simulated;
          Int.incr total_moves_simulated;
          let dst = apply_move m src in
          let cost = CubeTbl.find cost_so_far src + 1 in
          let dominated =
            match CubeTbl.find_opt cost_so_far dst with
            | Some c -> c <= cost
            | None -> false
          in
          if not dominated then begin
            CubeTbl.replace cost_so_far dst cost;
            CubeTbl.replace came_from dst (src, m);
            if is_goal dst then solution := Some (reconstruct came_from dst)
            else if !simulated < budget then
              let hw =
                if Float.(random_weight > 0.0) then
                  Float.max 0.01 (random_gauss 1.0 random_weight)
                else 1.0
              in
              let prio = Float.of_int cost +. (hw *. heuristic dst) in
              frontier := push_heap !frontier prio dst
          end
        end
      in

      while
        Option.is_none !solution
        && (not !frontier_exhausted)
        && !simulated < budget
      do
        match pop_heap !frontier with
        | None -> frontier_exhausted := true
        | Some (src, rest) ->
          frontier := rest;
          let last_move = Option.map ~f:snd (CubeTbl.find_opt came_from src) in
          for m = 0 to num_moves - 1 do
            if Option.is_none !solution && !simulated < budget then
              step_move src last_move m
          done
      done;

      match !solution with
      | Some s -> Some s
      | None when !frontier_exhausted || Float.(random_weight <= 0.0) -> None
      | None ->
        let tm = Unix.localtime (Unix.gettimeofday ()) in
        printf
          "[%02d:%02d:%02d] search budget of %d moves exceeded; restarting\n%!"
          tm.tm_hour tm.tm_min tm.tm_sec budget;
        attempt (Float.to_int (Float.of_int budget *. 1.5))
    in
    attempt max_moves

let count_solved pred (cube : cube) =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    if pred cubelets.(i) && is_cubelet_solved cubelets.(i) cube.(i) then
      Int.incr cnt
  done;
  !cnt

let count_bottom_edges_positioned (cube : cube) =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    let _, _, z = c in
    if z = -1 && norm1 c = 2 then
      if cube.(i) *@ (0, 0, -1) = (0, 0, -1) then Int.incr cnt
  done;
  !cnt

let count_bottom_corners_positioned (cube : cube) =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    let _, _, z = c in
    if z = -1 && norm1 c = 3 && is_cubelet_pos_solved c cube.(i) then
      Int.incr cnt
  done;
  !cnt

let log fmt =
  let tm = Unix.localtime (Unix.gettimeofday ()) in
  printf "[%02d:%02d:%02d] " tm.tm_hour tm.tm_min tm.tm_sec;
  printf (Stdlib.( ^^ ) fmt "\n%!")

let solve_layer name total is_goal heuristic rw cube =
  let curr = ref cube and moves_acc = ref [] in
  for i = 0 to total - 1 do
    log "%s #%d" name (i + 1);
    match astar !curr (is_goal i) (heuristic i) rw 100_000 with
    | Some (next_c, mvs) ->
      log "-> found solution with %d moves" (List.length mvs);
      curr := next_c;
      moves_acc := mvs :: !moves_acc
    | None -> failwith ("Failed " ^ name)
  done;
  (!curr, List.concat (List.rev !moves_acc))

let bottom_left_front_corner (cube : cube) =
  let target = (1, -1, -1) in
  let rec find i =
    if i = num_cubelets then failwith "Corner missing"
    else if cube.(i) *@ cubelets.(i) = target then (i, cube.(i))
    else find (i + 1)
  in
  find 0

let find_move n d =
  let rec loop i =
    if moves.(i).normal = n && moves.(i).dir = d then i else loop (i + 1)
  in
  loop 0

let solve_endgame (cube : cube) =
  let curr = ref cube and sol = ref [] in
  let left = find_move (0, -1, 0) 1 in
  let top = find_move (0, 0, 1) 1 in
  let bottom = find_move (0, 0, -1) 1 in
  let routine =
    [
      inv_move.(left);
      inv_move.(top);
      left;
      top;
      inv_move.(left);
      inv_move.(top);
      left;
      top;
    ]
  in
  let apply m =
    sol := m :: !sol;
    curr := apply_move m !curr
  in
  let is_corner_oriented () =
    let c_idx, r = bottom_left_front_corner !curr in
    let rot = ref r and solved = ref false in
    for _ = 0 to 3 do
      if is_cubelet_solved cubelets.(c_idx) !rot then solved := true;
      rot := rot_matrices.(bottom) *@* !rot
    done;
    !solved
  in
  for _ = 0 to 3 do
    while not (is_corner_oriented ()) do
      List.iter routine ~f:apply
    done;
    apply bottom
  done;
  while not (is_cube_solved !curr) do
    apply bottom
  done;
  (!curr, List.rev !sol)

let shuffle cube iters seed =
  Stdlib.Random.init seed;
  let curr = ref cube in
  for _ = 1 to iters do
    curr := apply_move (Stdlib.Random.int num_moves) !curr
  done;
  !curr

let solve (cube : cube) =
  let t0 = Unix.gettimeofday () in
  let start_sim = !total_moves_simulated in
  let c1, s1 =
    solve_layer "solving cubelet" 17
      (fun i c ->
        count_solved
          (fun v ->
            let _, _, z = v in
            z = 1 && norm1 v = 2
          )
          c
        >= Int.min 4 (i + 1)
        && count_solved (fun (_, _, z) -> z = 1) c >= Int.min 9 (i + 1)
        && count_solved (fun (_, _, z) -> z >= 0) c >= Int.min 17 (i + 1)
      )
      (fun i -> if i < 9 then top_layer_heuristic else middle_layer_heuristic)
      0.25 cube
  in
  log "--------------------------------------------------";
  let c2, s2 =
    solve_layer "solving bottom cross" 8
      (fun i c ->
        count_solved (fun (_, _, z) -> z >= 0) c = 17
        && count_bottom_edges_positioned c >= Int.min 4 (i + 1)
        && count_solved
             (fun v ->
               let _, _, z = v in
               z = -1 && norm1 v = 2
             )
             c
           >= Int.min 4 (i - 3)
      )
      (fun _ -> bottom_layer_edge_heuristic)
      0.25 c1
  in
  log "--------------------------------------------------";
  let c3, s3 =
    solve_layer "positioning bottom corners" 4
      (fun i c ->
        count_solved (fun (_, _, z) -> z >= 0) c = 17
        && count_solved
             (fun v ->
               let _, _, z = v in
               z = -1 && norm1 v = 2
             )
             c
           = 4
        && count_bottom_corners_positioned c >= Int.min 4 (i + 1)
      )
      (fun _ -> bottom_layer_corner_heuristic)
      0.30 c2
  in
  log "--------------------------------------------------";
  let c4, s4 = solve_endgame c3 in
  let moves = List.concat [ s1; s2; s3; s4 ] in
  let elapsed = Unix.gettimeofday () -. t0 in
  let moves_simulated = !total_moves_simulated - start_sim in
  log "Solved cube in %d moves." (List.length moves);
  log "is_cube_solved: %b" (is_cube_solved c4);
  log "- time elapsed: %.2f sec" elapsed;
  log "- moves simulated: %d (%.0f moves/sec)" moves_simulated
    (Float.of_int moves_simulated /. Float.max 0.001 elapsed);
  moves

let run_tests () =
  printf "Running OCaml Eigencube Invariant Tests...\n%!";
  (* Test 1: Cubelet count and canonical positions *)
  assert (num_cubelets = 26);
  assert (Array.length moves = 12);
  printf "  [PASS] Cubelet count and move counts match.\n%!";

  (* Test 2: Linear algebra and R * diag(c) = diag(c) invariant *)
  let c = (1, 1, 1) in
  assert (is_cubelet_solved c id3);
  let rot_x = rot_matrices.(0) in
  assert (not (is_cubelet_solved c rot_x));
  assert (dot (1, 0, 0) (0, 1, 0) = 0);
  assert (dot (1, 2, 3) (4, 5, 6) = 32);
  printf
    "  [PASS] Linear algebra and R * diag(c) = diag(c) invariant verified.\n%!";

  (* Test 3: 4x single-move identity *)
  let c0 = solved_cube () in
  for m = 0 to num_moves - 1 do
    let c1 = apply_move m c0 in
    let c2 = apply_move m c1 in
    let c3 = apply_move m c2 in
    let c4 = apply_move m c3 in
    assert (is_cube_solved c4);
    assert (not (is_cube_solved c1));
    assert (not (is_cube_solved c2));
    assert (not (is_cube_solved c3))
  done;
  printf "  [PASS] All 12 moves satisfy order-4 cyclic permutation.\n%!";

  (* Test 4: Inverse move cancellation *)
  for m = 0 to num_moves - 1 do
    let inv_m = inv_move.(m) in
    let c_after = apply_move inv_m (apply_move m c0) in
    assert (is_cube_solved c_after)
  done;
  printf "  [PASS] Inverse move cancellation verified for all 12 moves.\n%!";

  (* Test 5: 6x Sexy Move identity (R U R' U') * 6 = Identity *)
  let r = find_move (0, 1, 0) 1 in
  let u = find_move (0, 0, 1) 1 in
  let r_inv = inv_move.(r) in
  let u_inv = inv_move.(u) in
  let sexy = [ r; u; r_inv; u_inv ] in
  let c_ref = ref (solved_cube ()) in
  for _ = 1 to 6 do
    List.iter sexy ~f:(fun m -> c_ref := apply_move m !c_ref)
  done;
  assert (is_cube_solved !c_ref);
  printf "  [PASS] 6x Sexy Move ((R U R' U') * 6) restores identity.\n%!";

  (* Test 6: Deterministic search terminates with None when budget is exceeded *)
  let result = astar c0 (fun _ -> false) (fun _ -> 0.0) 0.0 500 in
  assert (Option.is_none result);
  printf
    "  [PASS] Budget exhaustion terminates cleanly with None under \
     deterministic search.\n\
     %!";

  (* Test 7: Short scramble solve round-trip *)
  let scrambled = shuffle (solved_cube ()) 10 999 in
  assert (not (is_cube_solved scrambled));
  let sol_moves = solve scrambled in
  let final_cube = ref scrambled in
  List.iter sol_moves ~f:(fun m -> final_cube := apply_move m !final_cube);
  assert (is_cube_solved !final_cube);
  printf "  [PASS] End-to-end solve verifies cube is completely solved.\n%!";

  (* Test 8: Successive solves isolation *)
  let c1 = shuffle (solved_cube ()) 5 123 in
  let c2 = shuffle (solved_cube ()) 5 456 in
  let m1 = solve c1 in
  let m2 = solve c2 in
  assert (List.length m1 > 0);
  assert (List.length m2 > 0);
  printf
    "  [PASS] Successive solves execute independently without interference.\n%!";

  printf "All OCaml invariant and solver tests passed successfully!\n%!"

let main () =
  let print_usage () =
    printf "Usage: dune exec ocaml/eigencube.exe -- [seed | --test]\n";
    printf
      "  [seed]   Solves a Rubik's cube scrambled from random seed (default: 42)\n";
    printf "  --test   Runs the complete invariant and solver test suite\n"
  in
  let args = Sys.get_argv () in
  if Array.length args > 1 then
    match args.(1) with
    | "--help" | "-h" -> print_usage ()
    | "--test" -> run_tests ()
    | s ->
    match Int.of_string_opt s with
    | Some seed ->
      log "Solving scrambled cube (seed=%d)..." seed;
      ignore (solve (shuffle (solved_cube ()) 100_000 seed))
    | None ->
      printf "Error: unrecognized option or invalid integer seed '%s'.\n\n%!" s;
      print_usage ();
      Stdlib.exit 1
  else begin
    log "Solving scrambled cube (seed=42)...";
    ignore (solve (shuffle (solved_cube ()) 100_000 42))
  end

let () = main ()
