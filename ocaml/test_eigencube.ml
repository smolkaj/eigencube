open Base
open Stdio
open Eigencube

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
  let c_sexy =
    Fn.apply_n_times ~n:6
      (fun c -> List.fold sexy ~init:c ~f:(fun acc m -> apply_move m acc))
      (solved_cube ())
  in
  assert (is_cube_solved c_sexy);
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
  let final_cube =
    List.fold sol_moves ~init:scrambled ~f:(fun acc m -> apply_move m acc)
  in
  assert (is_cube_solved final_cube);
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

let () = run_tests ()
