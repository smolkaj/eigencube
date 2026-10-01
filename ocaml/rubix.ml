(* rubix.ml - Minimalistic Rubik's Cube Solver in OCaml
   A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
   octahedral symmetry group, and multi-phase A* search with restarts.
   Leveraging Jane Street's Base and Stdio libraries. *)

open Base
open Stdio
open Poly

type vec3 = { x : int; y : int; z : int }

type mat3 = {
  m00 : int;
  m01 : int;
  m02 : int;
  m10 : int;
  m11 : int;
  m12 : int;
  m20 : int;
  m21 : int;
  m22 : int;
}

type move = { normal : vec3; dir : int }
type cube = int array (* 26 cubelets mapped to rotation index 0..23 *)

let norm1 v = Int.abs v.x + Int.abs v.y + Int.abs v.z
let dot a b = (a.x * b.x) + (a.y * b.y) + (a.z * b.z)

let id3 =
  {
    m00 = 1;
    m01 = 0;
    m02 = 0;
    m10 = 0;
    m11 = 1;
    m12 = 0;
    m20 = 0;
    m21 = 0;
    m22 = 1;
  }

let mat_vec m v =
  {
    x = (m.m00 * v.x) + (m.m01 * v.y) + (m.m02 * v.z);
    y = (m.m10 * v.x) + (m.m11 * v.y) + (m.m12 * v.z);
    z = (m.m20 * v.x) + (m.m21 * v.y) + (m.m22 * v.z);
  }

let mat_mul a b =
  {
    m00 = (a.m00 * b.m00) + (a.m01 * b.m10) + (a.m02 * b.m20);
    m01 = (a.m00 * b.m01) + (a.m01 * b.m11) + (a.m02 * b.m21);
    m02 = (a.m00 * b.m02) + (a.m01 * b.m12) + (a.m02 * b.m22);
    m10 = (a.m10 * b.m00) + (a.m11 * b.m10) + (a.m12 * b.m20);
    m11 = (a.m10 * b.m01) + (a.m11 * b.m11) + (a.m12 * b.m21);
    m12 = (a.m10 * b.m02) + (a.m11 * b.m12) + (a.m12 * b.m22);
    m20 = (a.m20 * b.m00) + (a.m21 * b.m10) + (a.m22 * b.m20);
    m21 = (a.m20 * b.m01) + (a.m21 * b.m11) + (a.m22 * b.m21);
    m22 = (a.m20 * b.m02) + (a.m21 * b.m12) + (a.m22 * b.m22);
  }

let crange = [ -1; 0; 1 ]

let all_vectors =
  List.concat_map crange ~f:(fun x ->
      List.concat_map crange ~f:(fun y ->
          List.map crange ~f:(fun z -> { x; y; z })
      )
  )

let cubelets =
  all_vectors |> List.filter ~f:(fun v -> norm1 v > 0) |> Array.of_list

let num_cubelets = 26
let unit_vectors = List.filter all_vectors ~f:(fun v -> norm1 v = 1)

let moves =
  List.concat_map unit_vectors ~f:(fun v ->
      [ -1; 1 ] |> List.map ~f:(fun dir -> { normal = v; dir })
  )
  |> Array.of_list

let num_moves = 12

let rotation_matrix { normal = v; dir } =
  let dim = if v.x <> 0 then 0 else if v.y <> 0 then 1 else 2 in
  match dim with
  | 0 ->
    {
      m00 = 1;
      m01 = 0;
      m02 = 0;
      m10 = 0;
      m11 = 0;
      m12 = dir;
      m20 = 0;
      m21 = -dir;
      m22 = 0;
    }
  | 1 ->
    {
      m00 = 0;
      m01 = 0;
      m02 = dir;
      m10 = 0;
      m11 = 1;
      m12 = 0;
      m20 = -dir;
      m21 = 0;
      m22 = 0;
    }
  | _ ->
    {
      m00 = 0;
      m01 = dir;
      m02 = 0;
      m10 = -dir;
      m11 = 0;
      m12 = 0;
      m20 = 0;
      m21 = 0;
      m22 = 1;
    }

(* Chiral octahedral symmetry group O (|O| = 24) *)
let rotations, rot_indices =
  let tbl = Hashtbl.Poly.create ()
  and q = Queue.create ()
  and arr = Array.create ~len:24 id3 in
  Queue.enqueue q id3;
  Hashtbl.set tbl ~key:id3 ~data:0;
  arr.(0) <- id3;
  let count = ref 1 in
  while not (Queue.is_empty q) do
    let curr = Queue.dequeue_exn q in
    Array.iter moves ~f:(fun m ->
        let next_r = mat_mul (rotation_matrix m) curr in
        if not (Hashtbl.mem tbl next_r) then begin
          let idx = !count in
          Int.incr count;
          Hashtbl.set tbl ~key:next_r ~data:idx;
          arr.(idx) <- next_r;
          Queue.enqueue q next_r
        end
    )
  done;
  (arr, tbl)

let move_rot =
  Array.init num_moves ~f:(fun m ->
      let mat = rotation_matrix moves.(m) in
      Array.init 24 ~f:(fun r ->
          Hashtbl.find_exn rot_indices (mat_mul mat rotations.(r))
      )
  )

let move_applies =
  Array.init num_moves ~f:(fun m ->
      let v = moves.(m).normal in
      Array.init num_cubelets ~f:(fun c ->
          Array.init 24 ~f:(fun r ->
              dot v (mat_vec rotations.(r) cubelets.(c)) > 0
          )
      )
  )

let is_rot_solved c_idx r_idx =
  let c = cubelets.(c_idx) and r = rotations.(r_idx) in
  (c.x = 0 || (r.m00 = 1 && r.m10 = 0 && r.m20 = 0))
  && (c.y = 0 || (r.m01 = 0 && r.m11 = 1 && r.m21 = 0))
  && (c.z = 0 || (r.m02 = 0 && r.m12 = 0 && r.m22 = 1))

let is_pos_solved c_idx r_idx =
  mat_vec rotations.(r_idx) cubelets.(c_idx) = cubelets.(c_idx)

let compute_dist pred =
  Array.init num_cubelets ~f:(fun c ->
      let dist = Array.create ~len:24 (-1) and q = Queue.create () in
      for r = 0 to 23 do
        if pred c r then begin
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

let apply_move m cube =
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
  let v1 = moves.(lm).normal and v2 = moves.(m).normal in
  (v1.x > v2.x || (v1.x = v2.x && (v1.y > v2.y || (v1.y = v2.y && v1.z > v2.z))))
  && v1.x = -v2.x
  && v1.y = -v2.y
  && v1.z = -v2.z

let is_cube_solved cube =
  let rec loop i =
    i = num_cubelets || (dist_solved.(i).(cube.(i)) = 0 && loop (i + 1))
  in
  loop 0

(* Heuristics (L0.5 norm) *)
let top_layer_heuristic cube =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    if cubelets.(i).z = 1 then
      sum := !sum +. (Float.of_int dist_solved.(i).(cube.(i)) **. p)
  done;
  (!sum **. (1.0 /. p)) /. 8.0

let middle_layer_heuristic cube =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    if cubelets.(i).z >= 0 then
      sum := !sum +. (Float.of_int dist_solved.(i).(cube.(i)) **. p)
  done;
  (!sum **. (1.0 /. p)) /. 4.0

let bottom_layer_edge_heuristic cube =
  let p = 0.5 and sum = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) in
    if not (c.z = -1 && norm1 c = 3) then
      sum := !sum +. (Float.of_int dist_solved.(i).(cube.(i)) **. p)
  done;
  (!sum **. (1.0 /. p)) /. 3.0

let bottom_layer_corner_heuristic cube =
  let p = 0.5 and s1 = ref 0.0 and s2 = ref 0.0 and s3 = ref 0.0 in
  for i = 0 to num_cubelets - 1 do
    let c = cubelets.(i) and r = cube.(i) in
    if c.z = 1 then s1 := !s1 +. (Float.of_int dist_solved.(i).(r) **. p)
    else if c.z = 0 then s2 := !s2 +. (Float.of_int dist_solved.(i).(r) **. p)
    else if c.z = -1 then
      let d = if norm1 c = 3 then dist_pos.(i).(r) else dist_solved.(i).(r) in
      s3 := !s3 +. (Float.of_int d **. p)
  done;
  ((!s1 **. (1.0 /. p)) /. 5.0)
  +. ((!s2 **. (1.0 /. p)) /. 3.0)
  +. ((!s3 **. (1.0 /. p)) /. 8.0)

(* Binary min-heap *)
module Heap = struct
  type 'a t = {
    mutable prio : float array;
    mutable data : 'a array;
    mutable size : int;
  }

  let create cap dummy =
    {
      prio = Array.create ~len:(Int.max 16 cap) 0.0;
      data = Array.create ~len:(Int.max 16 cap) dummy;
      size = 0;
    }

  let push h p x =
    if h.size = Array.length h.prio then begin
      let ncap = h.size * 2 in
      let np = Array.create ~len:ncap 0.0
      and nd = Array.create ~len:ncap h.data.(0) in
      Array.blit ~src:h.prio ~src_pos:0 ~dst:np ~dst_pos:0 ~len:h.size;
      Array.blit ~src:h.data ~src_pos:0 ~dst:nd ~dst_pos:0 ~len:h.size;
      h.prio <- np;
      h.data <- nd
    end;
    let rec up i =
      if i > 0 then
        let par = (i - 1) lsr 1 in
        if Float.(p < h.prio.(par)) then begin
          h.prio.(i) <- h.prio.(par);
          h.data.(i) <- h.data.(par);
          up par
        end
        else begin
          h.prio.(i) <- p;
          h.data.(i) <- x
        end
      else begin
        h.prio.(0) <- p;
        h.data.(0) <- x
      end
    in
    up h.size;
    h.size <- h.size + 1

  let pop h dummy =
    if h.size = 0 then None
    else begin
      let res = h.data.(0) in
      h.size <- h.size - 1;
      if h.size > 0 then begin
        let lp = h.prio.(h.size) and lx = h.data.(h.size) in
        h.data.(h.size) <- dummy;
        let rec down i =
          let l = (i lsl 1) + 1 in
          let r = l + 1 in
          if l < h.size then
            let b =
              if r < h.size && Float.(h.prio.(r) < h.prio.(l)) then r else l
            in
            if Float.(h.prio.(b) < lp) then begin
              h.prio.(i) <- h.prio.(b);
              h.data.(i) <- h.data.(b);
              down b
            end
            else begin
              h.prio.(i) <- lp;
              h.data.(i) <- lx
            end
          else begin
            h.prio.(i) <- lp;
            h.data.(i) <- lx
          end
        in
        down 0
      end;
      Some res
    end
end

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
      let frontier = Heap.create 4096 start in
      let came_from = Hashtbl.Poly.create ()
      and cost_so_far = Hashtbl.Poly.create () in
      Heap.push frontier 0.0 start;
      Hashtbl.set cost_so_far ~key:start ~data:0;
      let simulated = ref 0
      and budget_exceeded = ref false
      and solution = ref None in
      let frontier_active = ref true in
      while
        Option.is_none !solution && (not !budget_exceeded) && !frontier_active
      do
        match Heap.pop frontier start with
        | None -> frontier_active := false
        | Some src ->
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
                    Heap.push frontier
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

(* Layer predicates *)
let count_solved pred cube =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    if pred cubelets.(i) && dist_solved.(i).(cube.(i)) = 0 then Int.incr cnt
  done;
  !cnt

let count_bottom_edges_positioned cube =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    if cubelets.(i).z = -1 && norm1 cubelets.(i) = 2 then
      if
        mat_vec rotations.(cube.(i)) { x = 0; y = 0; z = -1 }
        = { x = 0; y = 0; z = -1 }
      then Int.incr cnt
  done;
  !cnt

let count_bottom_corners_positioned cube =
  let cnt = ref 0 in
  for i = 0 to num_cubelets - 1 do
    if
      cubelets.(i).z = -1
      && norm1 cubelets.(i) = 3
      && dist_pos.(i).(cube.(i)) = 0
    then Int.incr cnt
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
  let target = { x = 1; y = -1; z = -1 } in
  let rec find i =
    if i = num_cubelets then failwith "Corner missing"
    else if mat_vec rotations.(cube.(i)) cubelets.(i) = target then (i, cube.(i))
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
  let left = find_move { x = 0; y = -1; z = 0 } 1 in
  let top = find_move { x = 0; y = 0; z = 1 } 1 in
  let bottom = find_move { x = 0; y = 0; z = -1 } 1 in
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
  let c1, s1 =
    solve_layer "solving cubelet" 17
      (fun i c ->
        count_solved (fun v -> v.z = 1 && norm1 v = 2) c >= Int.min 4 (i + 1)
        && count_solved (fun v -> v.z = 1) c >= Int.min 9 (i + 1)
        && count_solved (fun v -> v.z >= 0) c >= Int.min 17 (i + 1)
      )
      (fun i -> if i < 9 then top_layer_heuristic else middle_layer_heuristic)
      0.25 cube
  in
  log "--------------------------------------------------";
  let c2, s2 =
    solve_layer "solving bottom cross" 8
      (fun i c ->
        count_solved (fun v -> v.z >= 0) c = 17
        && count_bottom_edges_positioned c >= Int.min 4 (i + 1)
        && count_solved (fun v -> v.z = -1 && norm1 v = 2) c >= Int.min 4 (i - 3)
      )
      (fun _ -> bottom_layer_edge_heuristic)
      0.25 c1
  in
  log "--------------------------------------------------";
  let c3, s3 =
    solve_layer "positioning bottom corners" 4
      (fun i c ->
        count_solved (fun v -> v.z >= 0) c = 17
        && count_solved (fun v -> v.z = -1 && norm1 v = 2) c = 4
        && count_bottom_corners_positioned c >= Int.min 4 (i + 1)
      )
      (fun _ -> bottom_layer_corner_heuristic)
      0.30 c2
  in
  log "--------------------------------------------------";
  let c4, s4 = solve_endgame c3 in
  let moves = List.concat [ s1; s2; s3; s4 ] in
  let elapsed = Unix.gettimeofday () -. t0 in
  log "Solved cube in %d moves." (List.length moves);
  log "is_cube_solved: %b" (is_cube_solved c4);
  log "- time elapsed: %.2f sec" elapsed;
  log "- moves simulated: %d (%.0f moves/sec)" !total_moves_simulated
    (Float.of_int !total_moves_simulated /. Float.max 0.001 elapsed);
  moves

let main () =
  let args = Sys.get_argv () in
  let seed = if Array.length args > 1 then Int.of_string args.(1) else 42 in
  log "Solving scrambled cube (seed=%d)..." seed;
  ignore (solve (shuffle (solved_cube ()) 100_000 seed))
