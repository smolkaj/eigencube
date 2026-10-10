// lib.rs - High-Performance Discrete Linear Algebra Rubik's Cube Solver (eigencube_perf)
// Exploits SO(3, Z) group structure (order 24), 26-byte states, L1-resident Cayley
// action tables, and precomputed heuristic lookups to achieve >8 million moves/sec.

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

// Compact 26-byte cube state: each byte is an index in 0..24 representing the cubelet's SO(3, Z) rotation
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct FastCube(pub [u8; 26]);

impl std::hash::Hash for FastCube {
  #[inline(always)]
  fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
    state.write(&self.0);
  }
}

pub const SOLVED_CUBE: FastCube = FastCube([0; 26]);

pub struct Tables {
  pub rotations: [Mat3; 24],
  pub next_rot: [[[u8; 24]; 26]; 12],
  pub sqrt_dist_solved: [[f32; 24]; 26],
  pub sqrt_target_dist: [[f32; 24]; 26],
  pub is_solved: [[bool; 24]; 26],
  pub is_bottom_edge_pos: [[bool; 24]; 26],
  pub is_bottom_corner_pos: [[bool; 24]; 26],
  pub prune_move: [[bool; 12]; 13],
  pub pos: [[Vec3; 24]; 26],
  pub gaussian_lut: [f32; 1024],
}

pub fn is_cubelet_solved(c: Vec3, r: Mat3) -> bool {
  (r * Mat3::diag(c)) == Mat3::diag(c)
}
pub fn is_cubelet_pos_solved(c: Vec3, r: Mat3) -> bool {
  (r * c) == c
}

impl Tables {
  pub fn get() -> &'static Tables {
    static INSTANCE: OnceLock<Tables> = OnceLock::new();
    INSTANCE.get_or_init(Self::init)
  }

  #[allow(clippy::needless_range_loop)]
  fn init() -> Self {
    let (mut rotations, mut rot_count, mut q, mut head, mut tail) = ([Mat3::ID; 24], 1, [0usize; 24], 0, 1);
    while head < tail {
      let r = rotations[q[head]];
      head += 1;
      for &m in &MOVES {
        let nxt = m.rot_mat() * r;
        if !rotations[..rot_count].contains(&nxt) {
          rotations[rot_count] = nxt;
          q[tail] = rot_count;
          tail += 1;
          rot_count += 1;
        }
      }
    }

    let mut next_rot = [[[0u8; 24]; 26]; 12];
    for m in 0..12 {
      let (mv, r_prime) = (MOVES[m], MOVES[m].rot_mat());
      for i in 0..26 {
        let c = CUBELETS[i];
        for r in 0..24 {
          let rot = rotations[r];
          next_rot[m][i][r] = if mv.normal.dot(rot * c) > 0 {
            rotations.iter().position(|&x| x == r_prime * rot).unwrap() as u8
          } else {
            r as u8
          };
        }
      }
    }

    let (mut is_solved, mut is_bottom_edge_pos, mut is_bottom_corner_pos, mut pos) =
      ([[false; 24]; 26], [[false; 24]; 26], [[false; 24]; 26], [[Vec3(0, 0, 0); 24]; 26]);
    for i in 0..26 {
      let c = CUBELETS[i];
      for r in 0..24 {
        let rot = rotations[r];
        is_solved[i][r] = is_cubelet_solved(c, rot);
        is_bottom_edge_pos[i][r] = c.is_bottom_edge() && (rot * Vec3(0, 0, -1)) == Vec3(0, 0, -1);
        is_bottom_corner_pos[i][r] = c.is_bottom_corner() && is_cubelet_pos_solved(c, rot);
        pos[i][r] = rot * c;
      }
    }

    let (mut sqrt_dist_solved, mut sqrt_target_dist) = ([[0.0f32; 24]; 26], [[0.0f32; 24]; 26]);
    for i in 0..26 {
      for is_target in [false, true] {
        if is_target && !CUBELETS[i].is_bottom_corner() {
          sqrt_target_dist[i] = sqrt_dist_solved[i];
          continue;
        }
        let is_goal = |r| if is_target { is_bottom_corner_pos[i][r] } else { is_solved[i][r] };
        let (mut d, mut queue, mut qh, mut qt) = ([usize::MAX; 24], [0usize; 24], 0, 0);
        for r in 0..24 {
          if is_goal(r) {
            d[r] = 0;
            queue[qt] = r;
            qt += 1;
          }
        }
        while qh < qt {
          let curr = queue[qh];
          qh += 1;
          for &m in &MOVES {
            let prev_mat = m.invert().rot_mat() * rotations[curr];
            let prev_r = rotations.iter().position(|&x| x == prev_mat).unwrap();
            if d[prev_r] == usize::MAX {
              d[prev_r] = d[curr] + 1;
              queue[qt] = prev_r;
              qt += 1;
            }
          }
        }
        let target = if is_target { &mut sqrt_target_dist[i] } else { &mut sqrt_dist_solved[i] };
        for r in 0..24 {
          target[r] = (d[r] as f32).sqrt();
        }
      }
    }

    let mut prune_move = [[false; 12]; 13];
    for last in 1..=12 {
      let p = MOVES[last - 1];
      for next in 0..12 {
        let m = MOVES[next];
        prune_move[last][next] =
          (p.normal == m.normal && p.dir == -m.dir) || (p.normal.dot(m.normal) == -1 && p.normal > m.normal);
      }
    }

    let mut gaussian_lut = [0.0f32; 1024];
    for k in 0..512 {
      let (u1, u2) = (((2 * k + 1) as f32) / 1024.0, (((k * 17 + 13) % 1024) as f32) / 1024.0);
      let (r, theta) = ((-2.0 * u1.ln()).sqrt(), 2.0 * std::f32::consts::PI * u2);
      gaussian_lut[2 * k] = r * theta.cos();
      gaussian_lut[2 * k + 1] = r * theta.sin();
    }

    Self {
      rotations,
      next_rot,
      sqrt_dist_solved,
      sqrt_target_dist,
      is_solved,
      is_bottom_edge_pos,
      is_bottom_corner_pos,
      prune_move,
      pos,
      gaussian_lut,
    }
  }
}

#[inline(always)]
pub fn apply_move(tables: &Tables, m: usize, cube: &FastCube) -> FastCube {
  let (mut next, table) = ([0u8; 26], &tables.next_rot[m]);
  for i in 0..26 {
    next[i] = table[i][cube.0[i] as usize];
  }
  FastCube(next)
}

pub fn is_cube_solved(tables: &Tables, cube: &FastCube) -> bool {
  (0..26).all(|i| tables.is_solved[i][cube.0[i] as usize])
}

const fn cubelets_where<const N: usize>(cond: u8) -> [usize; N] {
  let (mut arr, mut count, mut i) = ([0usize; N], 0, 0);
  while i < 26 {
    let c = CUBELETS[i];
    let ok = match cond {
      0 => c.2 == 1,
      1 => c.2 == 0,
      2 => c.2 >= 0,
      3 => c.2 == -1,
      4 => !c.is_bottom_corner(),
      5 => c.is_bottom_edge(),
      6 => c.is_bottom_corner(),
      _ => c.is_top_edge(),
    };
    if ok {
      arr[count] = i;
      count += 1;
    }
    i += 1;
  }
  arr
}

pub const TOP_CUBELETS: [usize; 9] = cubelets_where(0);
pub const MID_LAYER_CUBELETS: [usize; 8] = cubelets_where(1);
pub const MID_CUBELETS: [usize; 17] = cubelets_where(2);
pub const BOT_CUBELETS: [usize; 9] = cubelets_where(3);
pub const NON_BOTTOM_CORNER_CUBELETS: [usize; 22] = cubelets_where(4);
pub const BOTTOM_EDGE_CUBELETS: [usize; 4] = cubelets_where(5);
pub const BOTTOM_CORNER_CUBELETS: [usize; 4] = cubelets_where(6);
pub const TOP_EDGE_CUBELETS: [usize; 4] = cubelets_where(7);

#[inline(always)]
fn layer_heuristic(cube: &FastCube, indices: &[usize], table: &[[f32; 24]; 26], divisor: f32) -> f64 {
  let sum: f32 = indices.iter().map(|&i| table[i][cube.0[i] as usize]).sum();
  ((sum * sum) / divisor) as f64
}

#[inline(always)]
pub fn top_layer_heuristic(t: &Tables, c: &FastCube) -> f64 {
  layer_heuristic(c, &TOP_CUBELETS, &t.sqrt_dist_solved, 8.0)
}
#[inline(always)]
pub fn middle_layer_heuristic(t: &Tables, c: &FastCube) -> f64 {
  layer_heuristic(c, &MID_CUBELETS, &t.sqrt_dist_solved, 4.0)
}
#[inline(always)]
pub fn bottom_layer_edge_heuristic(t: &Tables, c: &FastCube) -> f64 {
  layer_heuristic(c, &NON_BOTTOM_CORNER_CUBELETS, &t.sqrt_dist_solved, 3.0)
}
#[inline(always)]
pub fn bottom_layer_corner_heuristic(t: &Tables, c: &FastCube) -> f64 {
  layer_heuristic(c, &TOP_CUBELETS, &t.sqrt_dist_solved, 5.0)
    + layer_heuristic(c, &MID_LAYER_CUBELETS, &t.sqrt_dist_solved, 3.0)
    + layer_heuristic(c, &BOT_CUBELETS, &t.sqrt_target_dist, 8.0)
}

pub fn log(msg: &str) {
  let s = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
  println!("[{:02}:{:02}:{:02}] {}", (s / 3600) % 24, (s / 60) % 60, s % 60, msg);
}

pub struct FastRng(pub u64);
impl FastRng {
  pub fn new(seed: u64) -> Self {
    Self(seed.max(1))
  }
  #[inline(always)]
  pub fn next_u16(&mut self) -> u16 {
    self.0 ^= self.0 << 13;
    self.0 ^= self.0 >> 7;
    self.0 ^= self.0 << 17;
    (self.0 >> 32) as u16
  }
}

// Fast A* Search with Flat 26-byte states and compact Priority Queue
pub fn astar(
  tables: &Tables,
  start: FastCube,
  is_goal: impl Fn(&FastCube) -> bool,
  heuristic: impl Fn(&FastCube) -> f64,
  random_weight: f64,
  max_moves: usize,
) -> Option<(FastCube, Vec<Move>)> {
  if is_goal(&start) {
    return Some((start, Vec::new()));
  }
  let (mut visited, mut rng, mut budget) = (FxHashMap::default(), FastRng::new(42), max_moves);

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
        if tables.prune_move[last_move as usize][m] {
          continue;
        }
        TOTAL_MOVES_SIMULATED.fetch_add(1, Ordering::Relaxed);
        simulated += 1;
        let (next_cost, dst) = (cost + 1, apply_move(tables, m, &state));
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
            path.push(MOVES[m_idx]);
            curr = apply_move(tables, m_idx ^ 1, &curr);
          }
          path.reverse();
          return Some((dst, path));
        }
        let hw = if random_weight > 0.0 {
          (1.0 + (random_weight as f32) * tables.gaussian_lut[(rng.next_u16() as usize) & 1023]).max(0.01) as f64
        } else {
          1.0
        };
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

fn count(cube: &FastCube, indices: &[usize], pred: impl Fn(usize, usize) -> bool) -> usize {
  indices.iter().filter(|&&i| pred(i, cube.0[i] as usize)).count()
}

fn top_done(tables: &Tables, i: usize, c: &FastCube) -> bool {
  count(c, &TOP_EDGE_CUBELETS, |i, r| tables.is_solved[i][r]) >= 4.min(i + 1)
    && count(c, &TOP_CUBELETS, |i, r| tables.is_solved[i][r]) >= 9.min(i + 1)
    && count(c, &MID_CUBELETS, |i, r| tables.is_solved[i][r]) >= 17.min(i + 1)
}
fn cross_done(tables: &Tables, i: usize, c: &FastCube) -> bool {
  top_done(tables, 16, c)
    && count(c, &BOTTOM_EDGE_CUBELETS, |i, r| tables.is_bottom_edge_pos[i][r]) >= 4.min(i + 1)
    && count(c, &BOTTOM_EDGE_CUBELETS, |i, r| tables.is_solved[i][r]) >= 4.min(i.saturating_sub(3))
}
fn corners_done(tables: &Tables, i: usize, c: &FastCube) -> bool {
  cross_done(tables, 7, c)
    && count(c, &BOTTOM_CORNER_CUBELETS, |i, r| tables.is_bottom_corner_pos[i][r]) >= 4.min(i + 1)
}

#[allow(clippy::too_many_arguments)]
fn solve_layer(
  tables: &Tables,
  name: &str,
  total: usize,
  is_goal: impl Fn(usize, &FastCube) -> bool,
  heuristic: impl Fn(usize, &FastCube) -> f64,
  random_weight: f64,
  cube: &mut FastCube,
  moves: &mut Vec<Move>,
) {
  for i in 0..total {
    log(&format!("{} #{}", name, i + 1));
    let (next_c, mvs) = astar(tables, *cube, |c| is_goal(i, c), |c| heuristic(i, c), random_weight, 100_000)
      .expect("layer search failed");
    log(&format!("-> found solution with {} moves", mvs.len()));
    *cube = next_c;
    moves.extend(mvs);
  }
  log("--------------------------------------------------");
}

// Endgame: corner orientation commutator and alignment
pub fn solve_endgame(tables: &Tables, cube: &mut FastCube, moves: &mut Vec<Move>) {
  let left_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, -1, 0) && m.dir == 1).unwrap();
  let top_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, 0, 1) && m.dir == 1).unwrap();
  let bot_idx = MOVES.iter().position(|&m| m.normal == Vec3(0, 0, -1) && m.dir == 1).unwrap();
  let routine = [left_idx ^ 1, top_idx ^ 1, left_idx, top_idx, left_idx ^ 1, top_idx ^ 1, left_idx, top_idx];
  let step = |c: &mut FastCube, mvs: &mut Vec<Move>, m: usize| {
    *c = apply_move(tables, m, c);
    mvs.push(MOVES[m]);
  };

  for _ in 0..4 {
    let corner_idx = (0..26).find(|&i| tables.pos[i][cube.0[i] as usize] == Vec3(1, -1, -1)).unwrap();
    while !(0..4).any(|rot_steps| {
      let mut r = cube.0[corner_idx] as usize;
      for _ in 0..rot_steps {
        r = tables.next_rot[bot_idx][corner_idx][r] as usize;
      }
      tables.is_solved[corner_idx][r]
    }) {
      for &m in &routine {
        step(cube, moves, m);
      }
    }
    step(cube, moves, bot_idx);
  }
  while !is_cube_solved(tables, cube) {
    step(cube, moves, bot_idx);
  }
}

pub fn shuffle(iters: usize, seed: u64, cube: FastCube) -> FastCube {
  let tables = Tables::get();
  let mut rng = crate::FastRng::new(seed);
  (0..iters).fold(cube, |c, _| apply_move(tables, (rng.next_u16() as usize) % 12, &c))
}

pub fn solve(mut cube: FastCube) -> Vec<Move> {
  let tables = Tables::get();
  let (t0, start_sim, mut moves) = (Instant::now(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed), Vec::new());
  let h1 = |i, c: &FastCube| if i < 9 { top_layer_heuristic(tables, c) } else { middle_layer_heuristic(tables, c) };
  let (h_cross, h_corners) = (
    |_, c: &FastCube| bottom_layer_edge_heuristic(tables, c),
    |_, c: &FastCube| bottom_layer_corner_heuristic(tables, c),
  );
  solve_layer(tables, "solving cubelet", 17, |i, c| top_done(tables, i, c), h1, 0.25, &mut cube, &mut moves);
  solve_layer(tables, "solving bottom cross", 8, |i, c| cross_done(tables, i, c), h_cross, 0.25, &mut cube, &mut moves);
  solve_layer(
    tables,
    "positioning bottom corners",
    4,
    |i, c| corners_done(tables, i, c),
    h_corners,
    0.30,
    &mut cube,
    &mut moves,
  );
  solve_endgame(tables, &mut cube, &mut moves);

  let (elapsed, sim) = (t0.elapsed().as_secs_f64(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed) - start_sim);
  log(&format!("Solved cube in {} moves.", moves.len()));
  log(&format!("is_cube_solved: {}", is_cube_solved(tables, &cube)));
  log(&format!("- time elapsed: {:.3} sec", elapsed));
  log(&format!("- moves simulated: {} ({:.0} moves/sec)", sim, sim as f64 / elapsed.max(0.0001)));
  moves
}
