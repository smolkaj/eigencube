use eigencube_perf::*;

#[test]
fn test_counts_and_group_size() {
  let tables = Tables::get();
  assert_eq!(CUBELETS.len(), 26);
  assert_eq!(MOVES.len(), 12);
  assert_eq!(tables.rotations.len(), 24);
  assert_eq!(tables.rotations[0], Mat3::ID);
}

#[test]
fn test_cubelet_index_groups() {
  for (i, &c) in CUBELETS.iter().enumerate() {
    assert_eq!(TOP_CUBELETS.contains(&i), c.2 == 1, "mismatch on TOP_CUBELETS for index {}", i);
    assert_eq!(MID_LAYER_CUBELETS.contains(&i), c.2 == 0, "mismatch on MID_LAYER_CUBELETS for index {}", i);
    assert_eq!(MID_CUBELETS.contains(&i), c.2 >= 0, "mismatch on MID_CUBELETS for index {}", i);
    assert_eq!(BOT_CUBELETS.contains(&i), c.2 == -1, "mismatch on BOT_CUBELETS for index {}", i);
    assert_eq!(TOP_EDGE_CUBELETS.contains(&i), c.is_top_edge(), "mismatch on TOP_EDGE_CUBELETS for index {}", i);
    assert_eq!(
      BOTTOM_EDGE_CUBELETS.contains(&i),
      c.is_bottom_edge(),
      "mismatch on BOTTOM_EDGE_CUBELETS for index {}",
      i
    );
    assert_eq!(
      BOTTOM_CORNER_CUBELETS.contains(&i),
      c.is_bottom_corner(),
      "mismatch on BOTTOM_CORNER_CUBELETS for index {}",
      i
    );
    assert_eq!(
      NON_BOTTOM_CORNER_CUBELETS.contains(&i),
      !c.is_bottom_corner(),
      "mismatch on NON_BOTTOM_CORNER_CUBELETS for index {}",
      i
    );
  }
}

#[test]
fn test_4x_single_move_identity() {
  let tables = Tables::get();
  let c0 = SOLVED_CUBE;
  for m in 0..12 {
    let c1 = apply_move(tables, m, &c0);
    let c2 = apply_move(tables, m, &c1);
    let c3 = apply_move(tables, m, &c2);
    let c4 = apply_move(tables, m, &c3);
    assert!(is_cube_solved(tables, &c4));
    assert!(!is_cube_solved(tables, &c1));
    assert!(!is_cube_solved(tables, &c2));
    assert!(!is_cube_solved(tables, &c3));
  }
}

#[test]
fn test_inverse_move_cancellation() {
  let tables = Tables::get();
  let c0 = SOLVED_CUBE;
  for m in 0..12 {
    let inv_m = m ^ 1;
    let c_after = apply_move(tables, inv_m, &apply_move(tables, m, &c0));
    assert!(is_cube_solved(tables, &c_after));
  }
}

#[test]
fn test_6x_sexy_move_identity() {
  let tables = Tables::get();
  let r_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, 1, 0) && m.dir == 1).unwrap();
  let u_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, 0, 1) && m.dir == 1).unwrap();
  let sexy = [r_idx, u_idx, r_idx ^ 1, u_idx ^ 1];

  let mut c = SOLVED_CUBE;
  for _ in 0..6 {
    for &m in &sexy {
      c = apply_move(tables, m, &c);
    }
  }
  assert!(is_cube_solved(tables, &c));
}

#[test]
fn test_deterministic_budget_exhaustion() {
  let tables = Tables::get();
  let c0 = SOLVED_CUBE;
  let result = astar(tables, c0, |_| false, |_| 0.0, 0.0, 10);
  assert!(result.is_none());
}

#[test]
fn test_full_scramble_solve() {
  let tables = Tables::get();
  let scrambled = shuffle(tables, 100, 42, SOLVED_CUBE);
  let moves = solve(scrambled);
  let mut curr = scrambled;
  for m in &moves {
    let m_idx = MOVES.iter().position(|x| x == m).unwrap();
    curr = apply_move(tables, m_idx, &curr);
  }
  assert!(is_cube_solved(tables, &curr));
}
