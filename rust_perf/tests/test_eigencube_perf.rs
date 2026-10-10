use eigencube_perf::*;

#[test]
fn test_counts_and_group_size() {
  let tables = Tables::get();
  assert_eq!(CUBELETS.len(), 26);
  assert_eq!(MOVES.len(), 12);
  assert_eq!(tables.rotations.len(), 24);
  assert_eq!(tables.rotations[0], Mat3::ID);
  for c_idx in 0..26 {
    for r_idx in 0..24 {
      assert!(tables.sqrt_dist_solved[c_idx][r_idx] <= 1.733);
      assert!(tables.sqrt_dist_pos[c_idx][r_idx] <= 1.733);
    }
  }
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
fn test_fast_cube_roundtrip() {
  let c0 = SOLVED_CUBE;
  let fast_cube = FastCube::from_cube(&c0);
  assert_eq!(fast_cube, FastCube::SOLVED);
  assert_eq!(fast_cube.to_cube(), c0);

  // Permuted cubelet ordering
  let mut permuted = c0;
  permuted.reverse();
  let fast_permuted = FastCube::from_cube(&permuted);
  assert_eq!(fast_permuted, FastCube::SOLVED);
  assert_eq!(fast_permuted.to_cube(), c0);

  let scrambled = shuffle(50, 42, c0);
  let fast_scrambled = FastCube::from_cube(&scrambled);
  assert_eq!(fast_scrambled.to_cube(), scrambled);
}

#[test]
fn test_4x_single_move_identity() {
  let c0 = FastCube::SOLVED;
  for m in 0..12 {
    let c1 = c0.apply_move(m);
    let c2 = c1.apply_move(m);
    let c3 = c2.apply_move(m);
    let c4 = c3.apply_move(m);
    assert!(c4.is_solved());
    assert!(!c1.is_solved());
    assert!(!c2.is_solved());
    assert!(!c3.is_solved());
  }
}

#[test]
fn test_inverse_move_cancellation() {
  let c0 = FastCube::SOLVED;
  for m in 0..12 {
    let inv_m = m ^ 1;
    let c_after = c0.apply_move(m).apply_move(inv_m);
    assert!(c_after.is_solved());
  }
}

#[test]
fn test_6x_sexy_move_identity() {
  let r_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, 1, 0) && m.dir == 1).unwrap();
  let u_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, 0, 1) && m.dir == 1).unwrap();
  let sexy = [r_idx, u_idx, r_idx ^ 1, u_idx ^ 1];

  let mut c = FastCube::SOLVED;
  for _ in 0..6 {
    for &m in &sexy {
      c = c.apply_move(m);
    }
  }
  assert!(c.is_solved());
}

#[test]
fn test_deterministic_budget_exhaustion() {
  let c0 = FastCube::SOLVED;
  let result = astar(c0, |_| false, |_| 0.0, 0.0, 10);
  assert!(result.is_none());
}

#[test]
fn test_full_scramble_solve() {
  let scrambled = shuffle(100, 42, SOLVED_CUBE);
  let moves = solve(scrambled);
  let mut curr = FastCube::from_cube(&scrambled);
  for m in &moves {
    let m_idx = MOVES.iter().position(|x| x == m).unwrap();
    curr = curr.apply_move(m_idx);
  }
  assert!(curr.is_solved());
}

#[test]
fn test_successive_solves_independent() {
  for seed in [1, 2] {
    let scrambled = shuffle(50, seed, SOLVED_CUBE);
    let moves = solve(scrambled);
    let mut curr = FastCube::from_cube(&scrambled);
    for m in &moves {
      let m_idx = MOVES.iter().position(|x| x == m).unwrap();
      curr = curr.apply_move(m_idx);
    }
    assert!(curr.is_solved());
  }
}
