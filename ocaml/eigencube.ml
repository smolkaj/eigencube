(* eigencube.ml - Minimalistic Rubik's Cube Solver in OCaml
   A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
   and multi-phase A* search with move-budgeted restarts. *)

open Base
open Stdio
open Poly

type vec = int * int * int [@@deriving compare, hash, sexp]
type mat = vec * vec * vec [@@deriving compare, hash, sexp]
type move = { normal : vec; dir : int } [@@deriving compare, sexp]

(* Hook required by [@@deriving hash] for array types in Base *)
let hash_fold_array f state arr = Array.fold arr ~init:state ~f

(* 26 cubelets, each mapped to its current 3x3 rotation matrix *)
module Cube = struct
  type t = mat array [@@deriving compare, hash, sexp]
end

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
      Array.findi_exn moves ~f:(fun _ mv -> mv = inv) |> fst
  )

(* Opposite face moves commute; prune duplicate branches by enforcing canonical order *)
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
let solved_cube () : Cube.t = Array.create ~len:num_cubelets id3

let apply_move m cube : Cube.t =
  let v = moves.(m).normal and rm = rot_matrices.(m) in
  Array.init num_cubelets ~f:(fun i ->
      let r = cube.(i) in
      if dot v (r *@ cubelets.(i)) > 0 then rm *@* r else r
  )

let is_cube_solved cube = Array.for_all2_exn cubelets cube ~f:is_cubelet_solved

let single_cubelet_bfs c r is_goal =
  if is_goal c r then 0
  else
    let q = Queue.create () in
    let visited = Hashtbl.Poly.create () in
    Queue.enqueue q (r, 0);
    Hashtbl.set visited ~key:r ~data:0;
    let rec bfs () =
      match Queue.dequeue q with
      | None -> 0
      | Some (curr_r, d) -> (
        let found =
          Array.find_map rot_matrices ~f:(fun rm ->
              let next_r = rm *@* curr_r in
              if Hashtbl.mem visited next_r then None
              else if is_goal c next_r then Some (d + 1)
              else begin
                Hashtbl.set visited ~key:next_r ~data:(d + 1);
                Queue.enqueue q (next_r, d + 1);
                None
              end
          )
        in
        match found with
        | Some ans -> ans
        | None -> bfs ()
      )
    in
    bfs ()

(* Lazy distance heuristics *)
let dist_solved_cache = Hashtbl.Poly.create ()

let min_moves_to_solved c r =
  Hashtbl.find_or_add dist_solved_cache (c, r) ~default:(fun () ->
      single_cubelet_bfs c r is_cubelet_solved
  )

let dist_pos_cache = Hashtbl.Poly.create ()

let min_moves_to_pos c r =
  Hashtbl.find_or_add dist_pos_cache (c, r) ~default:(fun () ->
      single_cubelet_bfs c r is_cubelet_pos_solved
  )

let top_layer_heuristic cube =
  let p = 0.5 in
  let sum =
    Array.fold2_exn cubelets cube ~init:0.0 ~f:(fun acc ((_, _, z) as c) r ->
        if z = 1 then acc +. (Float.of_int (min_moves_to_solved c r) **. p)
        else acc
    )
  in
  (sum **. (1.0 /. p)) /. 8.0

let middle_layer_heuristic cube =
  let p = 0.5 in
  let sum =
    Array.fold2_exn cubelets cube ~init:0.0 ~f:(fun acc ((_, _, z) as c) r ->
        if z >= 0 then acc +. (Float.of_int (min_moves_to_solved c r) **. p)
        else acc
    )
  in
  (sum **. (1.0 /. p)) /. 4.0

let bottom_layer_edge_heuristic cube =
  let p = 0.5 in
  let sum =
    Array.fold2_exn cubelets cube ~init:0.0 ~f:(fun acc ((_, _, z) as c) r ->
        if not (z = -1 && norm1 c = 3) then
          acc +. (Float.of_int (min_moves_to_solved c r) **. p)
        else acc
    )
  in
  (sum **. (1.0 /. p)) /. 3.0

let bottom_layer_corner_heuristic cube =
  let p = 0.5 in
  let s1, s2, s3 =
    Array.fold2_exn cubelets cube ~init:(0.0, 0.0, 0.0)
      ~f:(fun (s1, s2, s3) ((_, _, z) as c) r ->
        if z = 1 then
          (s1 +. (Float.of_int (min_moves_to_solved c r) **. p), s2, s3)
        else if z = 0 then
          (s1, s2 +. (Float.of_int (min_moves_to_solved c r) **. p), s3)
        else if z = -1 then
          let d =
            if norm1 c = 3 then min_moves_to_pos c r
            else min_moves_to_solved c r
          in
          (s1, s2, s3 +. (Float.of_int d **. p))
        else (s1, s2, s3)
    )
  in
  ((s1 **. (1.0 /. p)) /. 5.0)
  +. ((s2 **. (1.0 /. p)) /. 3.0)
  +. ((s3 **. (1.0 /. p)) /. 8.0)

(* Functional pairing heap priority queue *)
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
    match Hashtbl.find came_from curr with
    | Some (p, mv) -> loop p (mv :: acc)
    | None -> acc
  in
  (dst, loop dst [])

let should_prune last_move m =
  match last_move with
  | None -> false
  | Some lm -> m = inv_move.(lm) || opposite_pruned lm m

(* Multi-phase A* search with move-budgeted restarts (1.5x expansion) *)
let astar start is_goal heuristic random_weight max_moves =
  if is_goal start then Some (start, [])
  else
    let rec attempt budget =
      let frontier = ref (push_heap empty_heap 0.0 start) in
      let came_from = Hashtbl.create (module Cube) ~size:8192 in
      let cost_so_far = Hashtbl.create (module Cube) ~size:8192 in
      Hashtbl.set cost_so_far ~key:start ~data:0;
      let simulated = ref 0 in
      let frontier_exhausted = ref false in
      let solution = ref None in

      let step_move src last_move m =
        if not (should_prune last_move m) then begin
          Int.incr simulated;
          Int.incr total_moves_simulated;
          let dst = apply_move m src in
          let cost = Hashtbl.find_exn cost_so_far src + 1 in
          let dominated =
            match Hashtbl.find cost_so_far dst with
            | Some c -> c <= cost
            | None -> false
          in
          if not dominated then begin
            Hashtbl.set cost_so_far ~key:dst ~data:cost;
            Hashtbl.set came_from ~key:dst ~data:(src, m);
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
          let last_move = Option.map ~f:snd (Hashtbl.find came_from src) in
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

let count_solved pred cube =
  Array.counti cube ~f:(fun i r ->
      pred cubelets.(i) && is_cubelet_solved cubelets.(i) r
  )

let count_bottom_edges_positioned cube =
  Array.counti cube ~f:(fun i r ->
      let c = cubelets.(i) in
      let _, _, z = c in
      z = -1 && norm1 c = 2 && r *@ (0, 0, -1) = (0, 0, -1)
  )

let count_bottom_corners_positioned cube =
  Array.counti cube ~f:(fun i r ->
      let c = cubelets.(i) in
      let _, _, z = c in
      z = -1 && norm1 c = 3 && is_cubelet_pos_solved c r
  )

let log fmt =
  let tm = Unix.localtime (Unix.gettimeofday ()) in
  printf "[%02d:%02d:%02d] " tm.tm_hour tm.tm_min tm.tm_sec;
  printf (Stdlib.( ^^ ) fmt "\n%!")

let solve_layer name total is_goal heuristic rw cube =
  let rec loop i curr moves_acc =
    if i = total then (curr, List.concat (List.rev moves_acc))
    else begin
      log "%s #%d" name (i + 1);
      match astar curr (is_goal i) (heuristic i) rw 100_000 with
      | Some (next_c, mvs) ->
        log "-> found solution with %d moves" (List.length mvs);
        loop (i + 1) next_c (mvs :: moves_acc)
      | None -> failwith ("Failed " ^ name)
    end
  in
  loop 0 cube []

let bottom_left_front_corner cube =
  let target = (1, -1, -1) in
  Array.find_mapi_exn cube ~f:(fun i r ->
      if r *@ cubelets.(i) = target then Some (i, r) else None
  )

let find_move n d =
  Array.findi_exn moves ~f:(fun _ m -> m.normal = n && m.dir = d) |> fst

(* Endgame: orient bottom corners using (R' D' R D) * 2/4 and align bottom face *)
let solve_endgame cube =
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
  let apply (c, sol) m = (apply_move m c, m :: sol) in
  let apply_all state mvs = List.fold mvs ~init:state ~f:apply in
  let is_corner_oriented c =
    let c_idx, r = bottom_left_front_corner c in
    let bm = rot_matrices.(bottom) in
    let rec check k rot =
      k < 4
      && (is_cubelet_solved cubelets.(c_idx) rot || check (k + 1) (bm *@* rot))
    in
    check 0 r
  in
  let rec orient_corner state =
    if is_corner_oriented (fst state) then state
    else orient_corner (apply_all state routine)
  in
  let state =
    Fn.apply_n_times ~n:4 (fun s -> apply (orient_corner s) bottom) (cube, [])
  in
  let rec align_bottom state =
    if is_cube_solved (fst state) then state
    else align_bottom (apply state bottom)
  in
  let final_cube, rev_sol = align_bottom state in
  (final_cube, List.rev rev_sol)

let shuffle cube iters seed =
  Stdlib.Random.init seed;
  Fn.apply_n_times ~n:iters
    (fun c -> apply_move (Stdlib.Random.int num_moves) c)
    cube

(* Full 3-phase human solver: top layer -> middle edges -> bottom layer & endgame *)
let solve (cube : Cube.t) =
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
