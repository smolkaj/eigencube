(* rubix.ml - Minimalistic Rubik's Cube Solver in OCaml
   A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
   chiral octahedral symmetry group, and multi-phase A* search with restarts. *)

open Base
open Stdio
open Poly

type vec = int * int * int
type mat = vec * vec * vec
type move = { normal : vec; dir : int }
type cube = int array (* 26 cubelets mapped to rotation index 0..23 *)

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

let is_cubelet_solved c r =
  let colors = diag c in
  r *@* colors = colors

let is_cubelet_pos_solved c r = r *@ c = c

(* Chiral octahedral symmetry group O (|O| = 24) *)
let rotations, rot_indices =
  let tbl = Hashtbl.Poly.create () and q = Queue.create () in
  let arr = Array.create ~len:24 id3 in
  Queue.enqueue q id3;
  Hashtbl.set tbl ~key:id3 ~data:0;
  let count = ref 1 in
  while not (Queue.is_empty q) do
    let r = Queue.dequeue_exn q in
    Array.iter rot_matrices ~f:(fun rm ->
        let r' = rm *@* r in
        if not (Hashtbl.mem tbl r') then begin
          let idx = !count in
          Int.incr count;
          Hashtbl.set tbl ~key:r' ~data:idx;
          arr.(idx) <- r';
          Queue.enqueue q r'
        end
    )
  done;
  (arr, tbl)

let is_rot_solved c r = is_cubelet_solved cubelets.(c) rotations.(r)
let is_pos_solved c r = is_cubelet_pos_solved cubelets.(c) rotations.(r)

let move_rot =
  Array.init num_moves ~f:(fun m ->
      let rm = rot_matrices.(m) in
      Array.init 24 ~f:(fun r ->
          Hashtbl.find_exn rot_indices (rm *@* rotations.(r))
      )
  )

let move_applies =
  Array.init num_moves ~f:(fun m ->
      let v = moves.(m).normal in
      Array.init num_cubelets ~f:(fun c ->
          Array.init 24 ~f:(fun r -> dot v (rotations.(r) *@ cubelets.(c)) > 0)
      )
  )

let compute_dist is_goal =
  Array.init num_cubelets ~f:(fun c ->
      let dist = Array.create ~len:24 (-1) and q = Queue.create () in
      for r = 0 to 23 do
        if is_goal c r then begin
          dist.(r) <- 0;
          Queue.enqueue q r
        end
      done;
      while not (Queue.is_empty q) do
        let curr = Queue.dequeue_exn q in
        for m = 0 to num_moves - 1 do
          let nxt = move_rot.(m).(curr) in
          if dist.(nxt) = -1 then begin
            dist.(nxt) <- dist.(curr) + 1;
            Queue.enqueue q nxt
          end
        done
      done;
      dist
  )

let dist_solved = compute_dist is_rot_solved
let dist_pos = compute_dist is_pos_solved
let solved_cube () : cube = Array.create ~len:num_cubelets 0

let apply_move m (cube : cube) : cube =
  let res = Array.copy cube in
  let app = move_applies.(m) and tr = move_rot.(m) in
  for i = 0 to num_cubelets - 1 do
    let r = res.(i) in
    if app.(i).(r) then res.(i) <- tr.(r)
  done;
  res

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

let is_cube_solved (cube : cube) =
  let rec loop i =
    i = num_cubelets || (dist_solved.(i).(cube.(i)) = 0 && loop (i + 1))
  in
  loop 0

(* Heuristics (L0.5 norm) *)
let top_layer_heuristic cube =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let _, _, z = cubelets.(i) in
    if z = 1 then sum := !sum +. (Float.of_int dist_solved.(i).(cube.(i)) **. p)
  done;
  (!sum **. (1.0 /. p)) /. 8.0

let middle_layer_heuristic cube =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let _, _, z = cubelets.(i) in
    if z >= 0 then sum := !sum +. (Float.of_int dist_solved.(i).(cube.(i)) **. p)
  done;
  (!sum **. (1.0 /. p)) /. 4.0

let bottom_layer_edge_heuristic cube =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    let _, _, z = c in
    if not (z = -1 && norm1 c = 3) then
      sum := !sum +. (Float.of_int dist_solved.(i).(cube.(i)) **. p)
  done;
  (!sum **. (1.0 /. p)) /. 3.0

let bottom_layer_corner_heuristic cube =
  let p = 0.5 and s1 = ref 0.0 and s2 = ref 0.0 and s3 = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) and r = cube.(i) in
    let _, _, z = c in
    if z = 1 then s1 := !s1 +. (Float.of_int dist_solved.(i).(r) **. p)
    else if z = 0 then s2 := !s2 +. (Float.of_int dist_solved.(i).(r) **. p)
    else if z = -1 then
      let d = if norm1 c = 3 then dist_pos.(i).(r) else dist_solved.(i).(r) in
      s3 := !s3 +. (Float.of_int d **. p)
  done;
  ((!s1 **. (1.0 /. p)) /. 5.0)
  +. ((!s2 **. (1.0 /. p)) /. 3.0)
  +. ((!s3 **. (1.0 /. p)) /. 8.0)

(* Pairing heap for A* priority queue *)
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

let astar start is_goal heuristic random_weight max_moves =
  if is_goal start then Some (start, [])
  else
    let rec attempt budget =
      let frontier = ref (push_heap empty_heap 0.0 start) in
      let came_from = Hashtbl.Poly.create ()
      and cost_so_far = Hashtbl.Poly.create () in
      Hashtbl.set cost_so_far ~key:start ~data:0;
      let simulated = ref 0
      and budget_exceeded = ref false
      and active = ref true
      and solution = ref None in
      while Option.is_none !solution && (not !budget_exceeded) && !active do
        match pop_heap !frontier with
        | None -> active := false
        | Some (src, rest) ->
          frontier := rest;
          let last_m =
            match Hashtbl.find came_from src with
            | Some (_, m) -> Some m
            | None -> None
          in
          for m = 0 to num_moves - 1 do
            if Option.is_none !solution && not !budget_exceeded then begin
              let skip =
                match last_m with
                | Some lm -> m = inv_move.(lm) || opposite_pruned lm m
                | None -> false
              in
              if not skip then begin
                let dst = apply_move m src in
                Int.incr simulated;
                Int.incr total_moves_simulated;
                let cost = Hashtbl.find_exn cost_so_far src + 1 in
                if !simulated >= budget then budget_exceeded := true;
                match Hashtbl.find cost_so_far dst with
                | Some c when c <= cost -> ()
                | _ ->
                  Hashtbl.set cost_so_far ~key:dst ~data:cost;
                  Hashtbl.set came_from ~key:dst ~data:(src, m);
                  if is_goal dst then begin
                    let rec unwind curr acc =
                      match Hashtbl.find came_from curr with
                      | Some (p, mv) -> unwind p (mv :: acc)
                      | None -> acc
                    in
                    solution := Some (dst, unwind dst [])
                  end
                  else if not !budget_exceeded then begin
                    let hw =
                      if Float.(random_weight > 0.0) then
                        Float.max 0.01 (random_gauss 1.0 random_weight)
                      else 1.0
                    in
                    frontier :=
                      push_heap !frontier
                        (Float.of_int cost +. (hw *. heuristic dst))
                        dst
                  end
              end
            end
          done
      done;
      match !solution with
      | Some s -> Some s
      | None ->
        if (not !budget_exceeded) || Float.(random_weight <= 0.0) then None
        else begin
          let tm = Unix.localtime (Unix.gettimeofday ()) in
          printf
            "[%02d:%02d:%02d] search budget of %d moves exceeded; restarting\n\
             %!"
            tm.tm_hour tm.tm_min tm.tm_sec budget;
          attempt (Float.to_int (Float.of_int budget *. 1.5))
        end
    in
    attempt max_moves

let count_solved pred cube =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    if pred cubelets.(i) && dist_solved.(i).(cube.(i)) = 0 then Int.incr cnt
  done;
  !cnt

let count_bottom_edges_positioned cube =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    let _, _, z = c in
    if z = -1 && norm1 c = 2 then
      if rotations.(cube.(i)) *@ (0, 0, -1) = (0, 0, -1) then Int.incr cnt
  done;
  !cnt

let count_bottom_corners_positioned cube =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    let _, _, z = c in
    if z = -1 && norm1 c = 3 && dist_pos.(i).(cube.(i)) = 0 then Int.incr cnt
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
      moves_acc := List.append !moves_acc mvs
    | None -> failwith ("Failed " ^ name)
  done;
  (!curr, !moves_acc)

let bottom_left_front_corner cube =
  let target = (1, -1, -1) in
  let rec find i =
    if i = num_cubelets then failwith "Corner missing"
    else if rotations.(cube.(i)) *@ cubelets.(i) = target then (i, cube.(i))
    else find (i + 1)
  in
  find 0

let find_move n d =
  let rec loop i =
    if moves.(i).normal = n && moves.(i).dir = d then i else loop (i + 1)
  in
  loop 0

let solve_endgame cube =
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
    let c_idx, r_idx = bottom_left_front_corner !curr in
    let rot = ref r_idx and solved = ref false in
    for _ = 0 to 3 do
      if dist_solved.(c_idx).(!rot) = 0 then solved := true;
      rot := move_rot.(bottom).(!rot)
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

let solve cube =
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

let main () =
  let args = Sys.get_argv () in
  if
    Array.length args > 1
    && (String.equal args.(1) "--help" || String.equal args.(1) "-h")
  then begin
    printf "Usage: dune exec ocaml/main.exe -- [seed]\n";
    printf
      "Solves a Rubik's cube scrambled from the given random seed (default: 42).\n"
  end
  else begin
    let seed = if Array.length args > 1 then Int.of_string args.(1) else 42 in
    log "Solving scrambled cube (seed=%d)..." seed;
    ignore (solve (shuffle (solved_cube ()) 100_000 seed))
  end
