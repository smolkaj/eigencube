// lib.rs - High-Performance Discrete Linear Algebra Rubik's Cube Solver (eigencube_perf)
// Compiles SO(3, Z) transformations and distance heuristics into direct L1 lookup tables.

use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Mul;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;
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

pub const CUBELETS: [Vec3; 26] = {
  let mut cubelets = [Vec3(0, 0, 0); 26];
  let mut i = 0;
  while i < 26 {
    let idx = if i < 13 { i } else { i + 1 };
    cubelets[i] = Vec3((idx / 9) as i8 - 1, ((idx / 3) % 3) as i8 - 1, (idx % 3) as i8 - 1);
    i += 1;
  }
  cubelets
};

pub const UNIT_VECTORS: [Vec3; 6] =
  [Vec3(-1, 0, 0), Vec3(0, -1, 0), Vec3(0, 0, -1), Vec3(0, 0, 1), Vec3(0, 1, 0), Vec3(1, 0, 0)];

pub const MOVES: [Move; 12] = {
  let mut moves = [Move::new(Vec3(0, 0, 0), 0); 12];
  let mut i = 0;
  while i < 12 {
    moves[i] = Move::new(UNIT_VECTORS[i / 2], if i % 2 == 0 { -1 } else { 1 });
    i += 1;
  }
  moves
};

pub type Cubelet = (Vec3, Mat3);
pub type Cube = [Cubelet; 26];

pub const SOLVED_CUBE: Cube = {
  let mut cube = [(Vec3(0, 0, 0), Mat3::ID); 26];
  let mut i = 0;
  while i < 26 {
    cube[i] = (CUBELETS[i], Mat3::ID);
    i += 1;
  }
  cube
};

pub fn is_cubelet_solved(c: Vec3, r: Mat3) -> bool {
  (r * Mat3::diag(c)) == Mat3::diag(c)
}
pub fn is_cube_solved(cube: &Cube) -> bool {
  cube.iter().all(|&(c, r)| is_cubelet_solved(c, r))
}

const SQRT_TABLE: [f64; 5] = [0.0, 1.0, std::f64::consts::SQRT_2, 1.732_050_807_568_877_2, 2.0];

pub struct Tables {
  pub rotations: [Mat3; 24],
  pub transition: [[[u8; 24]; 26]; 12],
  pub solved_dist: [[u8; 24]; 26],
  pub pos_dist: [[u8; 24]; 26],
  pub prune_move: [[bool; 12]; 13],
}

impl Tables {
  pub fn get() -> &'static Tables {
    static INSTANCE: OnceLock<Tables> = OnceLock::new();
    INSTANCE.get_or_init(Self::init)
  }

  #[allow(clippy::needless_range_loop)]
  fn init() -> Self {
    let (mut rotations, mut count, mut i) = ([Mat3::ID; 24], 1, 0);
    while i < count {
      let r = rotations[i];
      for m in &MOVES {
        let nxt = m.rot_mat() * r;
        if !rotations[..count].contains(&nxt) {
          rotations[count] = nxt;
          count += 1;
        }
      }
      i += 1;
    }

    let mut transition = [[[0u8; 24]; 26]; 12];
    for m in 0..12 {
      let (mv, r_prime) = (MOVES[m], MOVES[m].rot_mat());
      for (c_idx, &c) in CUBELETS.iter().enumerate() {
        for (r_idx, &r) in rotations.iter().enumerate() {
          let next_r = if mv.normal.dot(r * c) > 0 { r_prime * r } else { r };
          transition[m][c_idx][r_idx] = rotations.iter().position(|&x| x == next_r).unwrap() as u8;
        }
      }
    }

    let compute_dist = |is_goal: fn(Vec3, Mat3) -> bool| {
      let mut dist = [[255u8; 24]; 26];
      for (c_idx, &c) in CUBELETS.iter().enumerate() {
        let (mut queue, mut qh, mut qt) = ([0usize; 24], 0, 0);
        for (r_idx, &r) in rotations.iter().enumerate() {
          if is_goal(c, r) {
            dist[c_idx][r_idx] = 0;
            queue[qt] = r_idx;
            qt += 1;
          }
        }
        while qh < qt {
          let curr = queue[qh];
          qh += 1;
          for m in 0..12 {
            let nxt_r = transition[m][c_idx][curr] as usize;
            if dist[c_idx][nxt_r] == 255 {
              dist[c_idx][nxt_r] = dist[c_idx][curr] + 1;
              queue[qt] = nxt_r;
              qt += 1;
            }
          }
        }
      }
      dist
    };

    let solved_dist = compute_dist(is_cubelet_solved);
    let pos_dist = compute_dist(|c, r| (r * c) == c);

    let mut prune_move = [[false; 12]; 13];
    for last in 1..=12 {
      for next in 0..12 {
        let (p, m) = (MOVES[last - 1], MOVES[next]);
        prune_move[last][next] =
          (p.normal == m.normal && p.dir == -m.dir) || (p.normal.dot(m.normal) == -1 && p.normal > m.normal);
      }
    }

    Self { rotations, transition, solved_dist, pos_dist, prune_move }
  }
}

// Compact 26-byte cube state: each byte is an index in 0..24 representing the cubelet's SO(3, Z) rotation
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct FastCube(pub [u8; 26]);

impl std::hash::Hash for FastCube {
  #[inline(always)]
  fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
    state.write(&self.0);
  }
}

impl FastCube {
  pub const SOLVED: Self = Self([0; 26]);

  pub fn from_cube(cube: &Cube) -> Self {
    let t = Tables::get();
    let mut arr = [0u8; 26];
    for (i, &(c, r)) in cube.iter().enumerate() {
      assert_eq!(c, CUBELETS[i]);
      arr[i] = t.rotations.iter().position(|&x| x == r).expect("invalid rotation") as u8;
    }
    Self(arr)
  }

  pub fn to_cube(&self) -> Cube {
    let t = Tables::get();
    let mut cube = SOLVED_CUBE;
    for i in 0..26 {
      cube[i] = (CUBELETS[i], t.rotations[self.0[i] as usize]);
    }
    cube
  }

  #[inline(always)]
  pub fn apply_move(&self, m: usize) -> Self {
    let trans = &Tables::get().transition[m];
    let mut next = [0u8; 26];
    for i in 0..26 {
      next[i] = trans[i][self.0[i] as usize];
    }
    Self(next)
  }

  #[inline(always)]
  pub fn is_solved(&self) -> bool {
    let t = Tables::get();
    (0..26).all(|i| t.solved_dist[i][self.0[i] as usize] == 0)
  }
}

#[inline(always)]
pub fn top_layer_heuristic(t: &Tables, c: &FastCube) -> f64 {
  let s: f64 = (0..26)
    .filter(|&i| CUBELETS[i].2 == 1)
    .map(|i| SQRT_TABLE[t.solved_dist[i][c.0[i] as usize] as usize])
    .sum();
  (s * s) / 8.0
}
#[inline(always)]
pub fn middle_layer_heuristic(t: &Tables, c: &FastCube) -> f64 {
  let s: f64 = (0..26)
    .filter(|&i| CUBELETS[i].2 >= 0)
    .map(|i| SQRT_TABLE[t.solved_dist[i][c.0[i] as usize] as usize])
    .sum();
  (s * s) / 4.0
}
#[inline(always)]
pub fn bottom_layer_edge_heuristic(t: &Tables, c: &FastCube) -> f64 {
  let s: f64 = (0..26)
    .filter(|&i| !CUBELETS[i].is_bottom_corner())
    .map(|i| SQRT_TABLE[t.solved_dist[i][c.0[i] as usize] as usize])
    .sum();
  (s * s) / 3.0
}
#[inline(always)]
pub fn bottom_layer_corner_heuristic(t: &Tables, c: &FastCube) -> f64 {
  let s_top: f64 = (0..26)
    .filter(|&i| CUBELETS[i].2 == 1)
    .map(|i| SQRT_TABLE[t.solved_dist[i][c.0[i] as usize] as usize])
    .sum();
  let s_mid: f64 = (0..26)
    .filter(|&i| CUBELETS[i].2 == 0)
    .map(|i| SQRT_TABLE[t.solved_dist[i][c.0[i] as usize] as usize])
    .sum();
  let s_bot: f64 = (0..26)
    .filter(|&i| CUBELETS[i].2 == -1)
    .map(|i| {
      let d =
        if CUBELETS[i].is_bottom_corner() { t.pos_dist[i][c.0[i] as usize] } else { t.solved_dist[i][c.0[i] as usize] };
      SQRT_TABLE[d as usize]
    })
    .sum();
  (s_top * s_top) / 5.0 + (s_mid * s_mid) / 3.0 + (s_bot * s_bot) / 8.0
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

// Fast A* Search with Flat 26-byte states, L1 lookup tables, and single-pass hash lookups
pub fn astar(
  start: FastCube,
  is_goal: impl Fn(&FastCube) -> bool,
  heuristic: impl Fn(&FastCube) -> f64,
  random_weight: f64,
  max_moves: usize,
) -> Option<(FastCube, Vec<usize>)> {
  if is_goal(&start) {
    return Some((start, Vec::new()));
  }
  let t = Tables::get();
  let (mut visited, mut rng, mut budget) =
    (FxHashMap::with_capacity_and_hasher(65536, Default::default()), Rng::new(42), max_moves);

  loop {
    visited.clear();
    visited.insert(start, (0u16, 0u8));
    let (mut frontier, mut simulated) = (BinaryHeap::from([(Reverse(0u64), 0u16, 0u8, start)]), 0);

    while simulated < budget {
      let Some((_, cost, last_move, state)) = frontier.pop() else { break };
      if cost > visited[&state].0 {
        continue;
      }
      for m in 0..12 {
        if simulated >= budget {
          break;
        }
        if t.prune_move[last_move as usize][m] {
          continue;
        }
        TOTAL_MOVES_SIMULATED.fetch_add(1, Ordering::Relaxed);
        simulated += 1;
        let (next_cost, dst) = (cost + 1, state.apply_move(m));
        let mv_code = (m + 1) as u8;

        match visited.entry(dst) {
          std::collections::hash_map::Entry::Occupied(mut e) => {
            if e.get().0 <= next_cost {
              continue;
            }
            e.insert((next_cost, mv_code));
          }
          std::collections::hash_map::Entry::Vacant(e) => {
            e.insert((next_cost, mv_code));
          }
        }

        if is_goal(&dst) {
          let (mut path, mut curr) = (Vec::new(), dst);
          while let Some(&(_, code)) = visited.get(&curr) {
            if code == 0 {
              break;
            }
            let m_idx = (code - 1) as usize;
            path.push(m_idx);
            curr = curr.apply_move(m_idx ^ 1);
          }
          path.reverse();
          return Some((dst, path));
        }

        let hw = if random_weight > 0.0 { rng.random_gauss(1.0, random_weight).max(0.01) } else { 1.0 };
        let prio = ((next_cost as f64) + hw * heuristic(&dst)).to_bits();
        frontier.push((Reverse(prio), next_cost, mv_code, dst));
      }
    }

    if simulated < budget || random_weight <= 0.0 {
      return None;
    }
    log(&format!("search budget of {} moves exceeded; restarting", budget));
    budget = (budget + 1).max((budget as f64 * 1.5) as usize);
  }
}

fn count(c: &FastCube, pred: impl Fn(usize) -> bool, is_solved: impl Fn(usize, usize) -> bool) -> usize {
  (0..26).filter(|&i| pred(i) && is_solved(i, c.0[i] as usize)).count()
}

fn top_done(t: &Tables, i: usize, c: &FastCube) -> bool {
  count(c, |j| CUBELETS[j].is_top_edge(), |j, r| t.solved_dist[j][r] == 0) >= 4.min(i + 1)
    && count(c, |j| CUBELETS[j].2 == 1, |j, r| t.solved_dist[j][r] == 0) >= 9.min(i + 1)
    && count(c, |j| CUBELETS[j].2 >= 0, |j, r| t.solved_dist[j][r] == 0) >= 17.min(i + 1)
}
fn cross_done(t: &Tables, i: usize, c: &FastCube) -> bool {
  top_done(t, 16, c)
    && count(c, |j| CUBELETS[j].is_bottom_edge(), |_, r| (t.rotations[r] * Vec3(0, 0, -1)) == Vec3(0, 0, -1))
      >= 4.min(i + 1)
    && count(c, |j| CUBELETS[j].is_bottom_edge(), |j, r| t.solved_dist[j][r] == 0) >= 4.min(i.saturating_sub(3))
}
fn corners_done(t: &Tables, i: usize, c: &FastCube) -> bool {
  cross_done(t, 7, c) && count(c, |j| CUBELETS[j].is_bottom_corner(), |j, r| t.pos_dist[j][r] == 0) >= 4.min(i + 1)
}

fn solve_layer(
  name: &str,
  total: usize,
  is_goal: impl Fn(usize, &FastCube) -> bool,
  heuristic: impl Fn(usize, &FastCube) -> f64,
  random_weight: f64,
  cube: &mut FastCube,
  moves: &mut Vec<usize>,
) {
  for i in 0..total {
    log(&format!("{} #{}", name, i + 1));
    let (next_c, mvs) =
      astar(*cube, |c| is_goal(i, c), |c| heuristic(i, c), random_weight, 100_000).expect("layer search failed");
    log(&format!("-> found solution with {} moves", mvs.len()));
    *cube = next_c;
    moves.extend(mvs);
  }
  log("--------------------------------------------------");
}

fn corner_oriented(t: &Tables, c_idx: usize, mut r_idx: usize, bot: usize) -> bool {
  (0..4).any(|_| {
    let ok = t.solved_dist[c_idx][r_idx] == 0;
    r_idx = t.transition[bot][c_idx][r_idx] as usize;
    ok
  })
}

pub fn solve_endgame(cube: &mut FastCube, moves: &mut Vec<usize>) {
  let t = Tables::get();
  let left = MOVES.iter().position(|&m| m.normal == Vec3(0, -1, 0) && m.dir == 1).unwrap();
  let top = MOVES.iter().position(|&m| m.normal == Vec3(0, 0, 1) && m.dir == 1).unwrap();
  let bot = MOVES.iter().position(|&m| m.normal == Vec3(0, 0, -1) && m.dir == 1).unwrap();
  let routine = [left ^ 1, top ^ 1, left, top, left ^ 1, top ^ 1, left, top];

  for _ in 0..4 {
    while !{
      let (c_idx, r_idx) = (0..26)
        .find_map(|i| {
          let r = cube.0[i] as usize;
          ((t.rotations[r] * CUBELETS[i]) == Vec3(1, -1, -1)).then_some((i, r))
        })
        .expect("corner not found");
      corner_oriented(t, c_idx, r_idx, bot)
    } {
      for &m in &routine {
        *cube = cube.apply_move(m);
        moves.push(m);
      }
    }
    *cube = cube.apply_move(bot);
    moves.push(bot);
  }
  while !cube.is_solved() {
    *cube = cube.apply_move(bot);
    moves.push(bot);
  }
}

pub fn shuffle(iters: usize, seed: u64, cube: Cube) -> Cube {
  let mut rng = Rng::new(seed);
  let mut fast_cube = FastCube::from_cube(&cube);
  for _ in 0..iters {
    let m = (rng.next_f64() * 12.0) as usize % 12;
    fast_cube = fast_cube.apply_move(m);
  }
  fast_cube.to_cube()
}

pub fn solve(cube: Cube) -> Vec<Move> {
  let (t0, start_sim) = (Instant::now(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed));
  let mut fast_cube = FastCube::from_cube(&cube);
  let mut moves = Vec::new();
  let t = Tables::get();

  let h1 = |i, c: &FastCube| if i < 9 { top_layer_heuristic(t, c) } else { middle_layer_heuristic(t, c) };
  let (h_cross, h_corners) =
    (|_, c: &FastCube| bottom_layer_edge_heuristic(t, c), |_, c: &FastCube| bottom_layer_corner_heuristic(t, c));

  solve_layer("solving cubelet", 17, |i, c| top_done(t, i, c), h1, 0.25, &mut fast_cube, &mut moves);
  solve_layer("solving bottom cross", 8, |i, c| cross_done(t, i, c), h_cross, 0.25, &mut fast_cube, &mut moves);
  solve_layer(
    "positioning bottom corners",
    4,
    |i, c| corners_done(t, i, c),
    h_corners,
    0.30,
    &mut fast_cube,
    &mut moves,
  );
  solve_endgame(&mut fast_cube, &mut moves);

  let (elapsed, sim) = (t0.elapsed().as_secs_f64(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed) - start_sim);
  log(&format!("Solved cube in {} moves.", moves.len()));
  log(&format!("is_cube_solved: {}", fast_cube.is_solved()));
  log(&format!("- time elapsed: {:.3} sec", elapsed));
  log(&format!("- moves simulated: {} ({:.0} moves/sec)", sim, sim as f64 / elapsed.max(0.0001)));
  moves.into_iter().map(|idx| MOVES[idx]).collect()
}
