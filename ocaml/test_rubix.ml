(* test_rubix.ml - Invariant and round-trip tests for the OCaml solver *)

open Base
open Stdio
open Poly

(* We can include or test core logic directly *)

(* Import rubix functions by compiling test_rubix *)
let () = printf "Running OCaml Rubix Invariant Tests...\n%!"

(* Test 1: Cubelet count and canonical positions *)
let test_cubelet_invariants () =
  assert (Rubix.num_cubelets = 26);
  assert (Array.length Rubix.moves = 12);
  printf "  [PASS] Cubelet count and move counts match.\n%!"

(* Test 2: 4x single-move identity *)
let test_four_turn_identity () =
  let c0 = Rubix.solved_cube () in
  for m = 0 to Rubix.num_moves - 1 do
    let c1 = Rubix.apply_move m c0 in
    let c2 = Rubix.apply_move m c1 in
    let c3 = Rubix.apply_move m c2 in
    let c4 = Rubix.apply_move m c3 in
    assert (Rubix.is_cube_solved c4);
    assert (not (Rubix.is_cube_solved c1));
    assert (not (Rubix.is_cube_solved c2));
    assert (not (Rubix.is_cube_solved c3))
  done;
  printf "  [PASS] All 12 moves satisfy order-4 cyclic permutation.\n%!"

(* Test 3: Inverse move cancellation *)
let test_inverse_cancellation () =
  let c0 = Rubix.solved_cube () in
  for m = 0 to Rubix.num_moves - 1 do
    let inv_m = Rubix.inv_move.(m) in
    let c_after = Rubix.apply_move inv_m (Rubix.apply_move m c0) in
    assert (Rubix.is_cube_solved c_after)
  done;
  printf "  [PASS] Inverse move cancellation verified for all 12 moves.\n%!"

(* Test 4: 6x Sexy Move identity (R U R' U') * 6 = Identity *)
let test_sexy_move_cycle () =
  let r = Rubix.find_move { x = 0; y = 1; z = 0 } 1 in
  let u = Rubix.find_move { x = 0; y = 0; z = 1 } 1 in
  let r_inv = Rubix.inv_move.(r) in
  let u_inv = Rubix.inv_move.(u) in
  let sexy = [ r; u; r_inv; u_inv ] in
  let c = ref (Rubix.solved_cube ()) in
  for _ = 1 to 6 do
    List.iter sexy ~f:(fun m -> c := Rubix.apply_move m !c)
  done;
  assert (Rubix.is_cube_solved !c);
  printf "  [PASS] 6x Sexy Move ((R U R' U') * 6) restores identity.\n%!"

(* Test 5: Unreachable goal terminates with None (no infinite restart loop) *)
let test_unreachable_goal () =
  let c = Rubix.solved_cube () in
  let result = Rubix.astar c (fun _ -> false) (fun _ -> 0.0) 0.0 500 in
  assert (Option.is_none result);
  printf
    "  [PASS] Unreachable goal terminates cleanly with None (no infinite \
     restart).\n\
     %!"

(* Test 6: Short scramble solve round-trip *)
let test_solve_roundtrip () =
  let scrambled = Rubix.shuffle (Rubix.solved_cube ()) 10 999 in
  assert (not (Rubix.is_cube_solved scrambled));
  let moves = Rubix.solve scrambled in
  let final_cube = ref scrambled in
  List.iter moves ~f:(fun m -> final_cube := Rubix.apply_move m !final_cube);
  assert (Rubix.is_cube_solved !final_cube);
  printf "  [PASS] End-to-end solve verifies cube is completely solved.\n%!"

let () =
  test_cubelet_invariants ();
  test_four_turn_identity ();
  test_inverse_cancellation ();
  test_sexy_move_cycle ();
  test_unreachable_goal ();
  test_solve_roundtrip ();
  printf "All OCaml invariant and solver tests passed successfully!\n%!"
