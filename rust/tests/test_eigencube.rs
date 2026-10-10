use eigencube::*;

#[test]
fn test_counts() {
  assert_eq!(CUBELETS.len(), 26);
  assert_eq!(MOVES.len(), 12);
}

#[test]
fn test_linear_algebra_invariants() {
  let c = Vec3(1, 1, 1);
  assert!(is_cubelet_solved(c, Mat3::ID));
  let rot_x = MOVES[0].rot_mat();
  assert!(!is_cubelet_solved(c, rot_x));
  assert_eq!(Vec3(1, 0, 0).dot(Vec3(0, 1, 0)), 0);
  assert_eq!(Vec3(1, 2, 3).dot(Vec3(4, 5, 6)), 32);
}

#[test]
fn test_4x_single_move_identity() {
  let c0 = SOLVED_CUBE;
  for &m in &MOVES {
    let c1 = apply_move(m, &c0);
    let c2 = apply_move(m, &c1);
    let c3 = apply_move(m, &c2);
    let c4 = apply_move(m, &c3);
    assert!(is_cube_solved(&c4));
    assert!(!is_cube_solved(&c1));
    assert!(!is_cube_solved(&c2));
    assert!(!is_cube_solved(&c3));
  }
}

#[test]
fn test_inverse_move_cancellation() {
  let c0 = SOLVED_CUBE;
  for &m in &MOVES {
    let inv_m = m.invert();
    let c_after = apply_move(inv_m, &apply_move(m, &c0));
    assert!(is_cube_solved(&c_after));
  }
}

#[test]
fn test_6x_sexy_move_identity() {
  let r = Move { normal: Vec3(0, 1, 0), dir: 1 };
  let u = Move { normal: Vec3(0, 0, 1), dir: 1 };
  let r_inv = r.invert();
  let u_inv = u.invert();
  let sexy = [r, u, r_inv, u_inv];

  let mut c = SOLVED_CUBE;
  for _ in 0..6 {
    for &m in &sexy {
      c = apply_move(m, &c);
    }
  }
  assert!(is_cube_solved(&c));
}

#[test]
fn test_deterministic_budget_exhaustion() {
  let c0 = SOLVED_CUBE;
  let result = astar(c0, |_| false, apply_move, |_| 0.0, 0.0, 10);
  assert!(result.is_none());
}

#[test]
fn test_unreachable_goal_frontier_exhaustion_with_random_weight() {
  // SO(3, Z) rotation group has 24 states. With budget = 100, frontier empties at 24 < budget.
  let result = astar(Mat3::ID, |_| false, |m, &r| m.rot_mat() * r, |_| 0.0, 0.25, 100);
  assert!(result.is_none());
}

#[test]
fn test_full_scramble_solve() {
  let scrambled = shuffle(100, 42, SOLVED_CUBE);
  let moves = solve(scrambled);
  let mut curr = scrambled;
  for &m in &moves {
    curr = apply_move(m, &curr);
  }
  assert!(is_cube_solved(&curr));
}

#[test]
fn test_successive_solves_independent() {
  for seed in [1, 2] {
    let scrambled = shuffle(50, seed, SOLVED_CUBE);
    let moves = solve(scrambled);
    let mut curr = scrambled;
    for &m in &moves {
      curr = apply_move(m, &curr);
    }
    assert!(is_cube_solved(&curr));
  }
}
