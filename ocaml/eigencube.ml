(* eigencube.ml - Minimalistic Rubik's Cube Solver in OCaml
   A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
   and multi-phase A* search with move-budgeted restarts. *)

open Base
open Stdio
open Poly

type vec = int * int * int [@@deriving compare, sexp]
type mat = vec * vec * vec [@@deriving compare, sexp]
type move = { normal : vec; dir : int } [@@deriving compare, sexp]
type cube = (vec * mat) list [@@deriving compare, sexp]

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

let all_vectors =
  let coords = [ -1; 0; 1 ] in
  List.Cartesian_product.map3 coords coords coords ~f:(fun x y z -> (x, y, z))

let cubelets = List.filter all_vectors ~f:(fun v -> norm1 v > 0)
let unit_vectors = List.filter all_vectors ~f:(fun v -> norm1 v = 1)

let moves =
  List.concat_map unit_vectors ~f:(fun normal ->
      [ { normal; dir = -1 }; { normal; dir = 1 } ]
  )

let rot_mat { normal = x, y, _; dir } =
  if x <> 0 then ((1, 0, 0), (0, 0, dir), (0, -dir, 0))
  else if y <> 0 then ((0, 0, dir), (0, 1, 0), (-dir, 0, 0))
  else ((0, dir, 0), (-dir, 0, 0), (0, 0, 1))

let invert_move { normal; dir } = { normal; dir = -dir }

let is_cubelet_solved c r =
  let colors = diag c in
  let sticker_directions = r *@* colors in
  sticker_directions = colors

let is_cube_solved (cube : cube) =
  List.for_all cube ~f:(fun (c, r) -> is_cubelet_solved c r)

let is_cubelet_pos_solved c r = r *@ c = c
let solved_cube : cube = List.map cubelets ~f:(fun c -> (c, id3))

let apply_move move (cube : cube) : cube =
  let v = move.normal in
  let r' = rot_mat move in
  List.map cube ~f:(fun ((c, r) as cubelet) ->
      if dot v (r *@ c) > 0 then (c, r' *@* r) else cubelet
  )

let total_moves_simulated = ref 0

let should_prune last_move move =
  match last_move with
  | None -> false
  | Some prev ->
    (prev.normal = move.normal && prev.dir = -move.dir)
    || (dot prev.normal move.normal = -1 && prev.normal > move.normal)

let log fmt =
  let tm = Unix.localtime (Unix.gettimeofday ()) in
  printf "[%02d:%02d:%02d] " tm.tm_hour tm.tm_min tm.tm_sec;
  printf (Stdlib.( ^^ ) fmt "\n%!")

(* Central Limit Theorem: sum of 12 uniform random floats has exact mean 6.0 and variance 1.0 *)
let random_gauss ~mean ~stdev =
  List.init 12 ~f:(fun _ -> Random.float 1.0)
  |> List.fold ~init:0.0 ~f:( +. )
  |> fun sum -> mean +. (stdev *. (sum -. 6.0))

(* Multi-phase A* search with move-budgeted restarts (1.5x expansion) *)
type 'a visited_node = { cost : int; prev : ('a * move) option }

type 'a frontier_node = {
  prio : float;
  cost : int;
  last_move : move option;
  state : 'a;
}

let reconstruct visited dst =
  let rec loop curr acc =
    match Hashtbl.find visited curr with
    | Some { prev = Some (p, mv); _ } -> loop p (mv :: acc)
    | _ -> acc
  in
  (dst, loop dst [])

let astar (type state) ~(start : state) ~(is_goal : state -> bool)
    ~(apply_move : move -> state -> state) ?(heuristic = fun _ -> 0.0)
    ?(random_weight = 0.0) ?(max_moves = 100_000) () =
  let visited = Hashtbl.Poly.create ~size:65536 () in
  let empty_frontier =
    Fheap.create ~compare:(fun a b -> Float.compare a.prio b.prio)
  in
  let rec attempt budget =
    Hashtbl.clear visited;
    Hashtbl.set visited ~key:start ~data:{ cost = 0; prev = None };
    search budget
      (Fheap.add empty_frontier
         { prio = 0.0; cost = 0; last_move = None; state = start }
      )
      0
  and search budget frontier simulated =
    match Fheap.pop frontier with
    | None -> None
    | Some ({ cost; last_move; state; _ }, rest_frontier) ->
    match
      expand_moves budget state (cost + 1) last_move moves rest_frontier
        simulated
    with
    | `Found solution -> Some solution
    | `Continue (next_frontier, next_simulated) ->
      if next_simulated < budget then search budget next_frontier next_simulated
      else if Float.(random_weight <= 0.0) then None
      else begin
        log "search budget of %d moves exceeded; restarting" budget;
        attempt
          (Int.max (budget + 1) (Float.to_int (Float.of_int budget *. 1.5)))
      end
  and expand_moves budget src next_cost last_move mvs frontier simulated =
    match mvs with
    | [] -> `Continue (frontier, simulated)
    | _ when simulated >= budget -> `Continue (frontier, simulated)
    | move :: rest when should_prune last_move move ->
      expand_moves budget src next_cost last_move rest frontier simulated
    | move :: rest -> (
      Int.incr total_moves_simulated;
      let dst = apply_move move src in
      match Hashtbl.find visited dst with
      | Some v when v.cost <= next_cost ->
        expand_moves budget src next_cost last_move rest frontier (simulated + 1)
      | _ ->
        Hashtbl.set visited ~key:dst
          ~data:{ cost = next_cost; prev = Some (src, move) };
        if is_goal dst then `Found (reconstruct visited dst)
        else
          let hw =
            if Float.(random_weight > 0.0) then
              Float.max 0.01 (random_gauss ~mean:1.0 ~stdev:random_weight)
            else 1.0
          in
          let prio = Float.of_int next_cost +. (hw *. heuristic dst) in
          expand_moves budget src next_cost last_move rest
            (Fheap.add frontier
               { prio; cost = next_cost; last_move = Some move; state = dst }
            )
            (simulated + 1)
    )
  in
  if is_goal start then Some (start, []) else attempt max_moves

(* Single-cubelet distance heuristics via unified A* search *)
let cubelet_dist ~cache ~is_goal c r =
  Hashtbl.find_or_add cache (c, r) ~default:(fun () ->
      match
        astar ~start:r ~is_goal:(is_goal c)
          ~apply_move:(fun m r -> rot_mat m *@* r)
          ()
      with
      | Some (_, path) -> List.length path
      | None -> 0
  )

let min_moves_to_solved =
  let cache = Hashtbl.Poly.create () in
  cubelet_dist ~cache ~is_goal:is_cubelet_solved

let min_moves_to_pos =
  let cache = Hashtbl.Poly.create () in
  cubelet_dist ~cache ~is_goal:is_cubelet_pos_solved

let norm_p05 cube ~cond ~f =
  List.sum
    (module Float)
    cube
    ~f:(fun (c, r) -> if cond c then Float.sqrt (Float.of_int (f c r)) else 0.0)
  |> fun sum -> sum *. sum

let top_layer_heuristic (cube : cube) =
  norm_p05 cube ~cond:(fun (_, _, z) -> z = 1) ~f:min_moves_to_solved /. 8.0

let middle_layer_heuristic (cube : cube) =
  norm_p05 cube ~cond:(fun (_, _, z) -> z >= 0) ~f:min_moves_to_solved /. 4.0

let bottom_layer_edge_heuristic (cube : cube) =
  norm_p05 cube
    ~cond:(fun ((_, _, z) as c) -> not (z = -1 && norm1 c = 3))
    ~f:min_moves_to_solved
  /. 3.0

let bottom_layer_corner_heuristic (cube : cube) =
  (norm_p05 cube ~cond:(fun (_, _, z) -> z = 1) ~f:min_moves_to_solved /. 5.0)
  +. (norm_p05 cube ~cond:(fun (_, _, z) -> z = 0) ~f:min_moves_to_solved /. 3.0)
  +. norm_p05 cube
       ~cond:(fun (_, _, z) -> z = -1)
       ~f:(fun c r ->
         if norm1 c = 3 then min_moves_to_pos c r else min_moves_to_solved c r
       )
     /. 8.0

let count_solved ~f (cube : cube) =
  List.count cube ~f:(fun (c, r) -> f c && is_cubelet_solved c r)

let count_bottom_edges_positioned (cube : cube) =
  List.count cube ~f:(fun (((_, _, z) as c), r) ->
      z = -1 && norm1 c = 2 && r *@ (0, 0, -1) = (0, 0, -1)
  )

let count_bottom_corners_positioned (cube : cube) =
  List.count cube ~f:(fun (((_, _, z) as c), r) ->
      z = -1 && norm1 c = 3 && is_cubelet_pos_solved c r
  )

let is_top_edge ((_, _, z) as v) = z = 1 && norm1 v = 2
let is_bottom_edge ((_, _, z) as v) = z = -1 && norm1 v = 2

let solve_layer ~name ~total ~is_goal ~heuristic ~random_weight cube =
  let rec loop i curr moves_acc =
    if i = total then (curr, List.concat (List.rev moves_acc))
    else begin
      log "%s #%d" name (i + 1);
      match
        astar ~start:curr ~is_goal:(is_goal i) ~apply_move
          ~heuristic:(heuristic i) ~random_weight ~max_moves:100_000 ()
      with
      | Some (next_c, mvs) ->
        log "-> found solution with %d moves" (List.length mvs);
        loop (i + 1) next_c (mvs :: moves_acc)
      | None -> failwith ("Failed " ^ name)
    end
  in
  loop 0 cube []

let bottom_left_front_corner (cube : cube) =
  List.find_exn cube ~f:(fun (c, r) -> r *@ c = (1, -1, -1))

(* Endgame: orient bottom corners using (R' D' R D) * 2/4 and align bottom face *)
let solve_endgame cube =
  let left = { normal = (0, -1, 0); dir = 1 } in
  let top = { normal = (0, 0, 1); dir = 1 } in
  let bottom = { normal = (0, 0, -1); dir = 1 } in
  let cycle = [ invert_move left; invert_move top; left; top ] in
  let routine = cycle @ cycle in
  let apply (c, sol) m = (apply_move m c, m :: sol) in
  let apply_all state mvs = List.fold mvs ~init:state ~f:apply in
  let is_corner_oriented c =
    let c_orig, r = bottom_left_front_corner c in
    let r' = rot_mat bottom in
    let rec check k r =
      k < 4 && (is_cubelet_solved c_orig r || check (k + 1) (r' *@* r))
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

let shuffle ~iters ~seed cube =
  Stdlib.Random.init seed;
  let moves_arr = Array.of_list moves in
  Fn.apply_n_times ~n:iters
    (fun c ->
      apply_move moves_arr.(Stdlib.Random.int (Array.length moves_arr)) c
    )
    cube

(* Full 3-phase human solver: top layer -> middle edges -> bottom layer & endgame *)
let solve (cube : cube) =
  let t0 = Unix.gettimeofday () in
  let start_sim = !total_moves_simulated in
  let c1, s1 =
    solve_layer ~name:"solving cubelet" ~total:17
      ~is_goal:(fun i c ->
        count_solved ~f:is_top_edge c >= Int.min 4 (i + 1)
        && count_solved ~f:(fun (_, _, z) -> z = 1) c >= Int.min 9 (i + 1)
        && count_solved ~f:(fun (_, _, z) -> z >= 0) c >= Int.min 17 (i + 1)
      )
      ~heuristic:(fun i ->
        if i < 9 then top_layer_heuristic else middle_layer_heuristic
      )
      ~random_weight:0.25 cube
  in
  log "--------------------------------------------------";
  let c2, s2 =
    solve_layer ~name:"solving bottom cross" ~total:8
      ~is_goal:(fun i c ->
        count_solved ~f:(fun (_, _, z) -> z >= 0) c = 17
        && count_bottom_edges_positioned c >= Int.min 4 (i + 1)
        && count_solved ~f:is_bottom_edge c >= Int.min 4 (i - 3)
      )
      ~heuristic:(fun _ -> bottom_layer_edge_heuristic)
      ~random_weight:0.25 c1
  in
  log "--------------------------------------------------";
  let c3, s3 =
    solve_layer ~name:"positioning bottom corners" ~total:4
      ~is_goal:(fun i c ->
        count_solved ~f:(fun (_, _, z) -> z >= 0) c = 17
        && count_solved ~f:is_bottom_edge c = 4
        && count_bottom_corners_positioned c >= Int.min 4 (i + 1)
      )
      ~heuristic:(fun _ -> bottom_layer_corner_heuristic)
      ~random_weight:0.30 c2
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
