// eigencube.rs - Minimalistic Rubik's Cube Solver in Rust
// A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
// and multi-phase A* search with move-budgeted restarts.

use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Mul;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

pub static TOTAL_MOVES_SIMULATED: AtomicUsize = AtomicUsize::new(0);

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Vec3(pub i8, pub i8, pub i8);
impl Vec3 {
  pub const fn norm1(self) -> i32 {
    self.0.abs() as i32 + self.1.abs() as i32 + self.2.abs() as i32
  }
  pub const fn dot(self, o: Vec3) -> i32 {
    (self.0 as i32 * o.0 as i32) + (self.1 as i32 * o.1 as i32) + (self.2 as i32 * o.2 as i32)
  }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Mat3(pub Vec3, pub Vec3, pub Vec3);
impl Mat3 {
  pub const ID: Mat3 = Mat3(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
  pub const fn diag(v: Vec3) -> Mat3 {
    Mat3(Vec3(v.0, 0, 0), Vec3(0, v.1, 0), Vec3(0, 0, v.2))
  }
  pub const fn transpose(self) -> Mat3 {
    let (a, b, c) = (self.0, self.1, self.2);
    Mat3(Vec3(a.0, b.0, c.0), Vec3(a.1, b.1, c.1), Vec3(a.2, b.2, c.2))
  }
}

impl Mul<Vec3> for Mat3 {
  type Output = Vec3;
  fn mul(self, v: Vec3) -> Vec3 {
    Vec3(self.0.dot(v) as i8, self.1.dot(v) as i8, self.2.dot(v) as i8)
  }
}

impl Mul<Mat3> for Mat3 {
  type Output = Mat3;
  fn mul(self, m: Mat3) -> Mat3 {
    let t = m.transpose();
    Mat3(self * t.0, self * t.1, self * t.2).transpose()
  }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Move {
  pub normal: Vec3,
  pub dir: i8,
}
impl Move {
  pub const fn invert(self) -> Self {
    Self { normal: self.normal, dir: -self.dir }
  }
  pub const fn rot_mat(self) -> Mat3 {
    let (Vec3(x, y, _), d) = (self.normal, self.dir);
    match (x != 0, y != 0) {
      (true, _) => Mat3(Vec3(1, 0, 0), Vec3(0, 0, d), Vec3(0, -d, 0)),
      (_, true) => Mat3(Vec3(0, 0, d), Vec3(0, 1, 0), Vec3(-d, 0, 0)),
      _ => Mat3(Vec3(0, d, 0), Vec3(-d, 0, 0), Vec3(0, 0, 1)),
    }
  }
}

pub type Cubelet = (Vec3, Mat3);
pub type Cube = [Cubelet; 26];

pub const CUBELETS: [Vec3; 26] = {
  let mut arr = [Vec3(0, 0, 0); 26];
  let mut i = 0;
  while i < 26 {
    let idx = if i < 13 { i } else { i + 1 };
    arr[i] = Vec3((idx / 9) as i8 - 1, ((idx / 3) % 3) as i8 - 1, (idx % 3) as i8 - 1);
    i += 1;
  }
  arr
};

pub const UNIT_VECTORS: [Vec3; 6] =
  [Vec3(-1, 0, 0), Vec3(0, -1, 0), Vec3(0, 0, -1), Vec3(0, 0, 1), Vec3(0, 1, 0), Vec3(1, 0, 0)];

pub const MOVES: [Move; 12] = {
  let mut arr = [Move { normal: Vec3(0, 0, 0), dir: 0 }; 12];
  let mut i = 0;
  while i < 12 {
    arr[i] = Move { normal: UNIT_VECTORS[i / 2], dir: if i % 2 == 0 { -1 } else { 1 } };
    i += 1;
  }
  arr
};

pub const SOLVED_CUBE: Cube = {
  let mut arr = [(Vec3(0, 0, 0), Mat3::ID); 26];
  let mut i = 0;
  while i < 26 {
    arr[i] = (CUBELETS[i], Mat3::ID);
    i += 1;
  }
  arr
};

pub fn is_cubelet_solved(c: Vec3, r: Mat3) -> bool {
  (r * Mat3::diag(c)) == Mat3::diag(c)
}
pub fn is_cube_solved(cube: &Cube) -> bool {
  cube.iter().all(|&(c, r)| is_cubelet_solved(c, r))
}
pub fn is_cubelet_pos_solved(c: Vec3, r: Mat3) -> bool {
  (r * c) == c
}

pub fn apply_move(m: Move, cube: &Cube) -> Cube {
  let r_prime = m.rot_mat();
  cube.map(|(c, r)| if m.normal.dot(r * c) > 0 { (c, r_prime * r) } else { (c, r) })
}

pub fn should_prune(last: Option<Move>, m: Move) -> bool {
  last
    .is_some_and(|p| (p.normal == m.normal && p.dir == -m.dir) || (p.normal.dot(m.normal) == -1 && p.normal > m.normal))
}

pub fn log(msg: &str) {
  let s = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
  println!("[{:02}:{:02}:{:02}] {}", (s / 3600) % 24, (s / 60) % 60, s % 60, msg);
}

pub struct Rng(pub u64);
impl Rng {
  pub fn new(seed: u64) -> Self {
    Self(seed.max(1))
  }
  pub fn next_f64(&mut self) -> f64 {
    self.0 ^= self.0 << 13;
    self.0 ^= self.0 >> 7;
    self.0 ^= self.0 << 17;
    (self.0 >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
  }
  pub fn random_gauss(&mut self, mean: f64, stdev: f64) -> f64 {
    mean + stdev * ((0..12).map(|_| self.next_f64()).sum::<f64>() - 6.0)
  }
}

// Multi-phase A* search with move-budgeted restarts (1.5x expansion)
pub fn astar<S: Copy + Ord + std::hash::Hash>(
  start: S,
  is_goal: impl Fn(&S) -> bool,
  apply_mv: impl Fn(Move, &S) -> S,
  heuristic: impl Fn(&S) -> f64,
  random_weight: f64,
  max_moves: usize,
) -> Option<(S, Vec<Move>)> {
  if is_goal(&start) {
    return Some((start, Vec::new()));
  }
  let mut visited: FxHashMap<S, (u32, Option<Move>)> = FxHashMap::default();
  visited.reserve(65536);
  let mut rng = Rng::new(42);
  let mut budget = max_moves;

  loop {
    visited.clear();
    visited.insert(start, (0, None));
    let mut frontier = BinaryHeap::new();
    frontier.push((Reverse(0u64), 0u32, None, start));
    let (mut simulated, mut goal_node) = (0, None);

    while goal_node.is_none() && simulated < budget {
      let Some((_, cost, last_move, state)) = frontier.pop() else { break };
      if cost > visited.get(&state).map_or(u32::MAX, |v| v.0) {
        continue;
      }
      for &m in &MOVES {
        if simulated >= budget || goal_node.is_some() || should_prune(last_move, m) {
          continue;
        }
        TOTAL_MOVES_SIMULATED.fetch_add(1, Ordering::Relaxed);
        simulated += 1;
        let next_cost = cost + 1;
        let dst = apply_mv(m, &state);

        if visited.get(&dst).is_some_and(|v| v.0 <= next_cost) {
          continue;
        }
        visited.insert(dst, (next_cost, Some(m)));
        if is_goal(&dst) {
          goal_node = Some(dst);
          break;
        }
        let hw = if random_weight > 0.0 { rng.random_gauss(1.0, random_weight).max(0.01) } else { 1.0 };
        let prio = (next_cost as f64) + hw * heuristic(&dst);
        frontier.push((Reverse(prio.to_bits()), next_cost, Some(m), dst));
      }
    }

    if let Some(target) = goal_node {
      let (mut path, mut curr) = (Vec::new(), target);
      while let Some(&(_, Some(m))) = visited.get(&curr) {
        path.push(m);
        curr = apply_mv(m.invert(), &curr);
      }
      path.reverse();
      return Some((target, path));
    }

    if random_weight <= 0.0 {
      return None;
    }
    log(&format!("search budget of {} moves exceeded; restarting", budget));
    budget = (budget + 1).max((budget as f64 * 1.5) as usize);
  }
}

type Cache = Mutex<Option<FxHashMap<(Vec3, Mat3), usize>>>;
static SOLVED_CACHE: Cache = Mutex::new(None);
static POS_CACHE: Cache = Mutex::new(None);

fn cubelet_dist(cache: &Cache, is_goal: fn(Vec3, Mat3) -> bool, c: Vec3, r: Mat3) -> usize {
  if let Some(&d) = cache.lock().unwrap().as_ref().and_then(|m| m.get(&(c, r))) {
    return d;
  }
  let d = astar(r, |&x| is_goal(c, x), |m, &x| m.rot_mat() * x, |_| 0.0, 0.0, 100_000).map_or(0, |(_, p)| p.len());
  cache.lock().unwrap().get_or_insert_with(FxHashMap::default).insert((c, r), d);
  d
}
pub fn min_moves_to_solved(c: Vec3, r: Mat3) -> usize {
  cubelet_dist(&SOLVED_CACHE, is_cubelet_solved, c, r)
}
pub fn min_moves_to_pos(c: Vec3, r: Mat3) -> usize {
  cubelet_dist(&POS_CACHE, is_cubelet_pos_solved, c, r)
}

pub fn norm_p05(cube: &Cube, cond: impl Fn(Vec3) -> bool, f: impl Fn(Vec3, Mat3) -> usize) -> f64 {
  let sum: f64 = cube.iter().filter(|&&(c, _)| cond(c)).map(|&(c, r)| (f(c, r) as f64).sqrt()).sum();
  sum * sum
}

pub fn top_layer_heuristic(cube: &Cube) -> f64 {
  norm_p05(cube, |c| c.2 == 1, min_moves_to_solved) / 8.0
}
pub fn middle_layer_heuristic(cube: &Cube) -> f64 {
  norm_p05(cube, |c| c.2 >= 0, min_moves_to_solved) / 4.0
}
pub fn bottom_layer_edge_heuristic(cube: &Cube) -> f64 {
  norm_p05(cube, |c| !is_bottom_corner(c), min_moves_to_solved) / 3.0
}
pub fn bottom_layer_corner_heuristic(cube: &Cube) -> f64 {
  (norm_p05(cube, |c| c.2 == 1, min_moves_to_solved) / 5.0)
    + (norm_p05(cube, |c| c.2 == 0, min_moves_to_solved) / 3.0)
    + (norm_p05(
      cube,
      |c| c.2 == -1,
      |c, r| if is_bottom_corner(c) { min_moves_to_pos(c, r) } else { min_moves_to_solved(c, r) },
    ) / 8.0)
}

pub fn count_solved(cube: &Cube, f: impl Fn(Vec3) -> bool) -> usize {
  cube.iter().filter(|&&(c, r)| f(c) && is_cubelet_solved(c, r)).count()
}
pub fn count_bottom_edges_positioned(cube: &Cube) -> usize {
  cube.iter().filter(|&&(c, r)| is_bottom_edge(c) && (r * Vec3(0, 0, -1)) == Vec3(0, 0, -1)).count()
}
pub fn count_bottom_corners_positioned(cube: &Cube) -> usize {
  cube.iter().filter(|&&(c, r)| is_bottom_corner(c) && is_cubelet_pos_solved(c, r)).count()
}
pub fn is_top_edge(v: Vec3) -> bool {
  v.2 == 1 && v.norm1() == 2
}
pub fn is_bottom_edge(v: Vec3) -> bool {
  v.2 == -1 && v.norm1() == 2
}
pub fn is_bottom_corner(v: Vec3) -> bool {
  v.2 == -1 && v.norm1() == 3
}

fn solve_layer(
  name: &str,
  total: usize,
  is_goal: impl Fn(usize, &Cube) -> bool,
  heuristic: impl Fn(usize) -> fn(&Cube) -> f64,
  random_weight: f64,
  mut cube: Cube,
  moves: &mut Vec<Move>,
) -> Cube {
  for i in 0..total {
    log(&format!("{} #{}", name, i + 1));
    let (next_c, mvs) =
      astar(cube, |c| is_goal(i, c), apply_move, heuristic(i), random_weight, 100_000).expect("layer search failed");
    log(&format!("-> found solution with {} moves", mvs.len()));
    cube = next_c;
    moves.extend(mvs);
  }
  log("--------------------------------------------------");
  cube
}

pub fn bottom_left_front_corner(cube: &Cube) -> Cubelet {
  *cube.iter().find(|&&(c, r)| (r * c) == Vec3(1, -1, -1)).expect("corner not found")
}

// Endgame: orient bottom corners using (R' D' R D) * 2/4 and align bottom face
pub fn solve_endgame(mut cube: Cube) -> (Cube, Vec<Move>) {
  let (l, u, d) = (
    Move { normal: Vec3(0, -1, 0), dir: 1 },
    Move { normal: Vec3(0, 0, 1), dir: 1 },
    Move { normal: Vec3(0, 0, -1), dir: 1 },
  );
  let cycle = [l.invert(), u.invert(), l, u];
  let routine = [cycle[0], cycle[1], cycle[2], cycle[3], cycle[0], cycle[1], cycle[2], cycle[3]];
  let mut solution = Vec::new();

  for _ in 0..4 {
    while !(0..4).any(|k| {
      let (c, r) = bottom_left_front_corner(&cube);
      is_cubelet_solved(c, (0..k).fold(r, |acc, _| d.rot_mat() * acc))
    }) {
      for &m in &routine {
        cube = apply_move(m, &cube);
        solution.push(m);
      }
    }
    cube = apply_move(d, &cube);
    solution.push(d);
  }
  while !is_cube_solved(&cube) {
    cube = apply_move(d, &cube);
    solution.push(d);
  }
  (cube, solution)
}

pub fn shuffle(iters: usize, seed: u64, mut cube: Cube) -> Cube {
  let mut rng = Rng::new(seed);
  for _ in 0..iters {
    let m = MOVES[(rng.next_f64() * MOVES.len() as f64) as usize % MOVES.len()];
    cube = apply_move(m, &cube);
  }
  cube
}

fn top_done(i: usize, c: &Cube) -> bool {
  count_solved(c, is_top_edge) >= 4.min(i + 1)
    && count_solved(c, |v| v.2 == 1) >= 9.min(i + 1)
    && count_solved(c, |v| v.2 >= 0) >= 17.min(i + 1)
}
fn cross_done(i: usize, c: &Cube) -> bool {
  top_done(16, c)
    && count_bottom_edges_positioned(c) >= 4.min(i + 1)
    && count_solved(c, is_bottom_edge) >= 4.min(i.saturating_sub(3))
}
fn corners_done(i: usize, c: &Cube) -> bool {
  cross_done(7, c) && count_bottom_corners_positioned(c) >= 4.min(i + 1)
}

// Full 3-phase human solver: top layer -> middle edges -> bottom layer & endgame
pub fn solve(cube: Cube) -> Vec<Move> {
  let t0 = Instant::now();
  let start_sim = TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed);
  let mut moves = Vec::new();

  let h1 = |i| if i < 9 { top_layer_heuristic } else { middle_layer_heuristic };
  let cube = solve_layer("solving cubelet", 17, top_done, h1, 0.25, cube, &mut moves);
  let cube =
    solve_layer("solving bottom cross", 8, cross_done, |_| bottom_layer_edge_heuristic, 0.25, cube, &mut moves);
  let cube = solve_layer(
    "positioning bottom corners",
    4,
    corners_done,
    |_| bottom_layer_corner_heuristic,
    0.30,
    cube,
    &mut moves,
  );

  let (c4, s4) = solve_endgame(cube);
  moves.extend(s4);

  let elapsed = t0.elapsed().as_secs_f64();
  let moves_simulated = TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed) - start_sim;
  log(&format!("Solved cube in {} moves.", moves.len()));
  log(&format!("is_cube_solved: {}", is_cube_solved(&c4)));
  log(&format!("- time elapsed: {:.2} sec", elapsed));
  log(&format!(
    "- moves simulated: {} ({:.0} moves/sec)",
    moves_simulated,
    moves_simulated as f64 / elapsed.max(0.001)
  ));
  moves
}
