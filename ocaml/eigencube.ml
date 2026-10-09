(* eigencube.ml - Minimalistic Rubik's Cube Solver in OCaml
   A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
   and multi-phase A* search with move-budgeted restarts. *)

open Base
open Stdio
open Poly

type vec = int * int * int [@@deriving compare, hash, sexp]
type mat = vec * vec * vec [@@deriving compare, hash, sexp]
type move = { normal : vec; dir : int } [@@deriving compare, sexp]

module Iarray = struct
  include Stdlib.Iarray

  let ( .%() ) = get
  let hash_fold_t f state arr = fold_left (fun s x -> f s x) state arr
  let sexp_of_t sexp_of_x arr = sexp_of_array sexp_of_x (to_array arr)
  let t_of_sexp x_of_sexp sexp = of_array (array_of_sexp x_of_sexp sexp)
end

let ( .%() ) = Iarray.get

(* 26 cubelets, each mapped to its current 3x3 rotation matrix *)
module Cube = struct
  type t = mat Iarray.t [@@deriving compare, hash, sexp]
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
  all_vectors |> List.filter ~f:(fun v -> norm1 v > 0) |> Iarray.of_list

let num_cubelets = Iarray.length cubelets
let unit_vectors = List.filter all_vectors ~f:(fun v -> norm1 v = 1)

let moves =
  unit_vectors
  |> List.concat_map ~f:(fun normal ->
      [ { normal; dir = -1 }; { normal; dir = 1 } ]
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

let is_cubelet_solved c r =
  let colors = diag c in
  let sticker_directions = r *@* colors in
  sticker_directions = colors

let is_cube_solved cube = Iarray.for_all2 is_cubelet_solved cubelets cube
let is_cubelet_pos_solved c r = r *@ c = c
let solved_cube : Cube.t = Iarray.init num_cubelets (fun _ -> id3)

let apply_move m cube : Cube.t =
  let v = moves.(m).normal and rm = rot_matrices.(m) in
  Iarray.init num_cubelets (fun i ->
      let r = cube.%(i) in
      if dot v (r *@ cubelets.%(i)) > 0 then rm *@* r else r
  )

let single_cubelet_bfs c r is_goal =
  With_return.with_return (fun { return } ->
      if is_goal c r then return 0;
      let q = Queue.create () in
      let visited = Hashtbl.Poly.create () in
      Queue.enqueue q (r, 0);
      Hashtbl.set visited ~key:r ~data:();
      while not (Queue.is_empty q) do
        let curr, d = Queue.dequeue_exn q in
        for i = 0 to num_moves - 1 do
          let next_r = rot_matrices.(i) *@* curr in
          if not (Hashtbl.mem visited next_r) then begin
            if is_goal c next_r then return (d + 1);
            Hashtbl.set visited ~key:next_r ~data:();
            Queue.enqueue q (next_r, d + 1)
          end
        done
      done;
      0
  )

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

let fold2 a b ~init ~f =
  let len = Iarray.length a in
  let rec loop i acc =
    if i = len then acc else loop (i + 1) (f acc a.%(i) b.%(i))
  in
  loop 0 init

let top_layer_heuristic cube =
  let p = 0.5 in
  let sum =
    fold2 cubelets cube ~init:0.0 ~f:(fun acc ((_, _, z) as c) r ->
        if z = 1 then acc +. (Float.of_int (min_moves_to_solved c r) **. p)
        else acc
    )
  in
  (sum **. (1.0 /. p)) /. 8.0

let middle_layer_heuristic cube =
  let p = 0.5 in
  let sum =
    fold2 cubelets cube ~init:0.0 ~f:(fun acc ((_, _, z) as c) r ->
        if z >= 0 then acc +. (Float.of_int (min_moves_to_solved c r) **. p)
        else acc
    )
  in
  (sum **. (1.0 /. p)) /. 4.0

let bottom_layer_edge_heuristic cube =
  let p = 0.5 in
  let sum =
    fold2 cubelets cube ~init:0.0 ~f:(fun acc ((_, _, z) as c) r ->
        if not (z = -1 && norm1 c = 3) then
          acc +. (Float.of_int (min_moves_to_solved c r) **. p)
        else acc
    )
  in
  (sum **. (1.0 /. p)) /. 3.0

let bottom_layer_corner_heuristic cube =
  let p = 0.5 in
  let s1, s2, s3 =
    fold2 cubelets cube ~init:(0.0, 0.0, 0.0)
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

let empty_frontier =
  Fheap.create ~compare:(fun (p1, _) (p2, _) -> Float.compare p1 p2)

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

(* Prune branches that immediately invert the previous move, or that apply
   commuting moves on opposite faces (e.g. L R vs R L) in non-canonical order. *)
let should_prune last_move move =
  match last_move with
  | None -> false
  | Some prev ->
    move = inv_move.(prev)
    ||
    let n1 = moves.(prev).normal and n2 = moves.(move).normal in
    dot n1 n2 = -1 && Poly.(n1 > n2)

(* Multi-phase A* search with move-budgeted restarts (1.5x expansion) *)
let astar start is_goal heuristic random_weight max_moves =
  let rec attempt budget =
    let came_from = Hashtbl.create (module Cube) ~size:8192 in
    let cost_so_far = Hashtbl.create (module Cube) ~size:8192 in
    Hashtbl.set cost_so_far ~key:start ~data:0;

    let rec expand_moves src last_move move frontier simulated =
      if move = num_moves || simulated >= budget then
        `Continue (frontier, simulated)
      else if should_prune last_move move then
        expand_moves src last_move (move + 1) frontier simulated
      else begin
        Int.incr total_moves_simulated;
        let simulated = simulated + 1 in
        let dst = apply_move move src in
        let cost = Hashtbl.find_exn cost_so_far src + 1 in
        let dominated =
          match Hashtbl.find cost_so_far dst with
          | Some c -> c <= cost
          | None -> false
        in
        if dominated then
          expand_moves src last_move (move + 1) frontier simulated
        else begin
          Hashtbl.set cost_so_far ~key:dst ~data:cost;
          Hashtbl.set came_from ~key:dst ~data:(src, move);
          if is_goal dst then `Found (reconstruct came_from dst)
          else
            let hw =
              if Float.(random_weight > 0.0) then
                Float.max 0.01 (random_gauss 1.0 random_weight)
              else 1.0
            in
            let prio = Float.of_int cost +. (hw *. heuristic dst) in
            expand_moves src last_move (move + 1)
              (Fheap.add frontier (prio, dst))
              simulated
        end
      end
    in

    let rec search frontier simulated =
      match Fheap.pop frontier with
      | None -> None
      | Some ((_prio, src), rest_frontier) -> (
        let last_move = Option.map ~f:snd (Hashtbl.find came_from src) in
        match expand_moves src last_move 0 rest_frontier simulated with
        | `Found solution -> Some solution
        | `Continue (next_frontier, next_simulated) ->
          if next_simulated >= budget then
            if Float.(random_weight <= 0.0) then None
            else begin
              let tm = Unix.localtime (Unix.gettimeofday ()) in
              printf
                "[%02d:%02d:%02d] search budget of %d moves exceeded; restarting\n\
                 %!"
                tm.tm_hour tm.tm_min tm.tm_sec budget;
              attempt (Float.to_int (Float.of_int budget *. 1.5))
            end
          else search next_frontier next_simulated
      )
    in
    search (Fheap.add empty_frontier (0.0, start)) 0
  in
  if is_goal start then Some (start, []) else attempt max_moves

let count_matching f cube =
  let len = Iarray.length cube in
  let rec loop i acc =
    if i = len then acc
    else
      let inc = if f cubelets.%(i) cube.%(i) then 1 else 0 in
      loop (i + 1) (acc + inc)
  in
  loop 0 0

let count_solved pred cube =
  count_matching (fun c r -> pred c && is_cubelet_solved c r) cube

let count_bottom_edges_positioned cube =
  count_matching
    (fun c r ->
      let _, _, z = c in
      z = -1 && norm1 c = 2 && r *@ (0, 0, -1) = (0, 0, -1)
    )
    cube

let count_bottom_corners_positioned cube =
  count_matching
    (fun c r ->
      let _, _, z = c in
      z = -1 && norm1 c = 3 && is_cubelet_pos_solved c r
    )
    cube

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
  let len = Iarray.length cube in
  let rec loop i =
    if i = len then failwith "Corner not found"
    else
      let r = cube.%(i) in
      if r *@ cubelets.%(i) = target then (i, r) else loop (i + 1)
  in
  loop 0

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
      && (is_cubelet_solved cubelets.%(c_idx) rot || check (k + 1) (bm *@* rot))
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
