// lib.rs - Minimalistic Rubik's Cube Solver in Rust (eigencube crate)
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
  pub const fn is_top_edge(self) -> bool {
    self.2 == 1 && self.norm1() == 2
  }
  pub const fn is_bottom_edge(self) -> bool {
    self.2 == -1 && self.norm1() == 2
  }
  pub const fn is_bottom_corner(self) -> bool {
    self.2 == -1 && self.norm1() == 3
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
  pub const fn new(normal: Vec3, dir: i8) -> Self {
    Self { normal, dir }
  }
  pub const fn invert(self) -> Self {
    Self::new(self.normal, -self.dir)
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
pub const UNIT_VECTORS: [Vec3; 6] =
  [Vec3(-1, 0, 0), Vec3(0, -1, 0), Vec3(0, 0, -1), Vec3(0, 0, 1), Vec3(0, 1, 0), Vec3(1, 0, 0)];
const fn init_tables() -> ([Vec3; 26], Cube, [Move; 12]) {
  let (mut c, mut s, mut m, mut i) =
    ([Vec3(0, 0, 0); 26], [(Vec3(0, 0, 0), Mat3::ID); 26], [Move::new(Vec3(0, 0, 0), 0); 12], 0);
  while i < 26 {
    let idx = if i < 13 { i } else { i + 1 };
    c[i] = Vec3((idx / 9) as i8 - 1, ((idx / 3) % 3) as i8 - 1, (idx % 3) as i8 - 1);
    s[i] = (c[i], Mat3::ID);
    if i < 12 {
      m[i] = Move::new(UNIT_VECTORS[i / 2], if i % 2 == 0 { -1 } else { 1 });
    }
    i += 1;
  }
  (c, s, m)
}
pub const CUBELETS: [Vec3; 26] = init_tables().0;
pub const SOLVED_CUBE: Cube = init_tables().1;
pub const MOVES: [Move; 12] = init_tables().2;

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
  let Some(p) = last else { return false };
  (p.normal == m.normal && p.dir == -m.dir) || (p.normal.dot(m.normal) == -1 && p.normal > m.normal)
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
  let (mut visited, mut rng, mut budget) = (FxHashMap::default(), Rng::new(42), max_moves);

  loop {
    visited.clear();
    visited.insert(start, (0, None));
    let (mut frontier, mut simulated) = (BinaryHeap::from([(Reverse(0), 0, None, start)]), 0);

    while simulated < budget {
      let Some((_, cost, last_move, state)) = frontier.pop() else { break };
      if cost > visited[&state].0 {
        continue;
      }
      for &m in &MOVES {
        if simulated >= budget {
          break;
        }
        if should_prune(last_move, m) {
          continue;
        }
        TOTAL_MOVES_SIMULATED.fetch_add(1, Ordering::Relaxed);
        simulated += 1;
        let (next_cost, dst) = (cost + 1, apply_mv(m, &state));

        if visited.get(&dst).is_some_and(|v| v.0 <= next_cost) {
          continue;
        }
        visited.insert(dst, (next_cost, Some(m)));
        if is_goal(&dst) {
          let (mut path, mut curr) = (Vec::new(), dst);
          while let Some(&(_, Some(mv))) = visited.get(&curr) {
            path.push(mv);
            curr = apply_mv(mv.invert(), &curr);
          }
          path.reverse();
          return Some((dst, path));
        }
        let hw = if random_weight > 0.0 { rng.random_gauss(1.0, random_weight).max(0.01) } else { 1.0 };
        frontier.push((Reverse(((next_cost as f64) + hw * heuristic(&dst)).to_bits()), next_cost, Some(m), dst));
      }
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
  cube.iter().filter(|&&(c, _)| cond(c)).map(|&(c, r)| (f(c, r) as f64).sqrt()).sum::<f64>().powi(2)
}

pub fn top_layer_heuristic(cube: &Cube) -> f64 {
  norm_p05(cube, |c| c.2 == 1, min_moves_to_solved) / 8.0
}
pub fn middle_layer_heuristic(cube: &Cube) -> f64 {
  norm_p05(cube, |c| c.2 >= 0, min_moves_to_solved) / 4.0
}
pub fn bottom_layer_edge_heuristic(cube: &Cube) -> f64 {
  norm_p05(cube, |c| !c.is_bottom_corner(), min_moves_to_solved) / 3.0
}
pub fn bottom_layer_corner_heuristic(cube: &Cube) -> f64 {
  let c_h = |c: Vec3, r| if c.is_bottom_corner() { min_moves_to_pos(c, r) } else { min_moves_to_solved(c, r) };
  (norm_p05(cube, |c| c.2 == 1, min_moves_to_solved) / 5.0)
    + (norm_p05(cube, |c| c.2 == 0, min_moves_to_solved) / 3.0)
    + (norm_p05(cube, |c| c.2 == -1, c_h) / 8.0)
}

pub fn count_solved(cube: &Cube, f: impl Fn(Vec3) -> bool) -> usize {
  cube.iter().filter(|&&(c, r)| f(c) && is_cubelet_solved(c, r)).count()
}
pub fn count_bottom_edges_positioned(cube: &Cube) -> usize {
  cube.iter().filter(|&&(c, r)| c.is_bottom_edge() && (r * Vec3(0, 0, -1)) == Vec3(0, 0, -1)).count()
}
pub fn count_bottom_corners_positioned(cube: &Cube) -> usize {
  cube.iter().filter(|&&(c, r)| c.is_bottom_corner() && is_cubelet_pos_solved(c, r)).count()
}

fn solve_layer(
  name: &str,
  total: usize,
  is_goal: impl Fn(usize, &Cube) -> bool,
  heuristic: impl Fn(usize, &Cube) -> f64,
  random_weight: f64,
  cube: &mut Cube,
  moves: &mut Vec<Move>,
) {
  for i in 0..total {
    log(&format!("{} #{}", name, i + 1));
    let (next_c, mvs) = astar(*cube, |c| is_goal(i, c), apply_move, |c| heuristic(i, c), random_weight, 100_000)
      .expect("layer search failed");
    log(&format!("-> found solution with {} moves", mvs.len()));
    *cube = next_c;
    moves.extend(mvs);
  }
  log("--------------------------------------------------");
}

pub fn bottom_left_front_corner(cube: &Cube) -> Cubelet {
  *cube.iter().find(|&&(c, r)| (r * c) == Vec3(1, -1, -1)).expect("corner not found")
}

fn step(cube: &mut Cube, moves: &mut Vec<Move>, m: Move) {
  *cube = apply_move(m, cube);
  moves.push(m);
}
fn corner_oriented((c, mut r): Cubelet, d_rot: Mat3) -> bool {
  (0..4).any(|_| {
    let ok = is_cubelet_solved(c, r);
    r = d_rot * r;
    ok
  })
}
// Endgame: orient bottom corners using (R' D' R D) * 2/4 and align bottom face
pub fn solve_endgame(cube: &mut Cube, moves: &mut Vec<Move>) {
  let (l, u, d) = (Move::new(Vec3(0, -1, 0), 1), Move::new(Vec3(0, 0, 1), 1), Move::new(Vec3(0, 0, -1), 1));
  let (routine, d_rot) = ([l.invert(), u.invert(), l, u, l.invert(), u.invert(), l, u], d.rot_mat());

  for _ in 0..4 {
    while !corner_oriented(bottom_left_front_corner(cube), d_rot) {
      for &m in &routine {
        step(cube, moves, m);
      }
    }
    step(cube, moves, d);
  }
  while !is_cube_solved(cube) {
    step(cube, moves, d);
  }
}

pub fn shuffle(iters: usize, seed: u64, cube: Cube) -> Cube {
  let mut rng = Rng::new(seed);
  (0..iters).fold(cube, |c, _| apply_move(MOVES[(rng.next_f64() * 12.0) as usize % 12], &c))
}

fn top_done(i: usize, c: &Cube) -> bool {
  count_solved(c, Vec3::is_top_edge) >= 4.min(i + 1)
    && count_solved(c, |v| v.2 == 1) >= 9.min(i + 1)
    && count_solved(c, |v| v.2 >= 0) >= 17.min(i + 1)
}
fn cross_done(i: usize, c: &Cube) -> bool {
  top_done(16, c)
    && count_bottom_edges_positioned(c) >= 4.min(i + 1)
    && count_solved(c, Vec3::is_bottom_edge) >= 4.min(i.saturating_sub(3))
}
fn corners_done(i: usize, c: &Cube) -> bool {
  cross_done(7, c) && count_bottom_corners_positioned(c) >= 4.min(i + 1)
}

// Full 3-phase human solver: top layer -> middle edges -> bottom layer & endgame
pub fn solve(mut cube: Cube) -> Vec<Move> {
  let (t0, start_sim, mut moves) = (Instant::now(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed), Vec::new());

  let h1 = |i, c: &Cube| if i < 9 { top_layer_heuristic(c) } else { middle_layer_heuristic(c) };
  let (h_cross, h_corners) =
    (|_, c: &Cube| bottom_layer_edge_heuristic(c), |_, c: &Cube| bottom_layer_corner_heuristic(c));
  solve_layer("solving cubelet", 17, top_done, h1, 0.25, &mut cube, &mut moves);
  solve_layer("solving bottom cross", 8, cross_done, h_cross, 0.25, &mut cube, &mut moves);
  solve_layer("positioning bottom corners", 4, corners_done, h_corners, 0.30, &mut cube, &mut moves);
  solve_endgame(&mut cube, &mut moves);

  let (elapsed, sim) = (t0.elapsed().as_secs_f64(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed) - start_sim);
  log(&format!("Solved cube in {} moves.", moves.len()));
  log(&format!("is_cube_solved: {}", is_cube_solved(&cube)));
  log(&format!("- time elapsed: {:.2} sec", elapsed));
  log(&format!("- moves simulated: {} ({:.0} moves/sec)", sim, sim as f64 / elapsed.max(0.001)));
  moves
}
