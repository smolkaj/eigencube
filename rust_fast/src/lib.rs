use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::{hash_map::Entry, BinaryHeap};
use std::ops::Mul;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;
use std::time::Instant;

pub static TOTAL_MOVES_SIMULATED: AtomicUsize = AtomicUsize::new(0);

// -- Discrete 3D Linear Algebra: SO(3, ℤ) rotations & Euclidean vectors -------

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Vec3(pub i8, pub i8, pub i8);

impl Vec3 {
  pub const fn dot(self, o: Vec3) -> i32 {
    self.0 as i32 * o.0 as i32 + self.1 as i32 * o.1 as i32 + self.2 as i32 * o.2 as i32
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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
  fn mul(self, o: Mat3) -> Mat3 {
    let t = o.transpose();
    Mat3(self * t.0, self * t.1, self * t.2).transpose()
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Move {
  pub normal: Vec3,
  pub dir: i8,
}

impl Move {
  /// Exact closed-form Rodrigues rotation formula for 90-degree turns of unit axis v:
  /// R = v vᵀ - d [v]ₓ, where [v]ₓ is the cross-product skew-symmetric matrix.
  /// https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  pub const fn rot_mat(self) -> Mat3 {
    let (Vec3(x, y, z), d) = (self.normal, self.dir);
    Mat3(Vec3(x * x, d * z, -d * y), Vec3(-d * z, y * y, d * x), Vec3(d * y, -d * x, z * z))
  }
}

pub type Cubelet = (Vec3, Mat3);
pub type Cube = [Cubelet; 26];

pub const CUBELETS: [Vec3; 26] = {
  let (mut c, mut i) = ([Vec3(0, 0, 0); 26], 0);
  while i < 26 {
    let j = if i < 13 { i } else { i + 1 };
    c[i] = Vec3((j / 9) as i8 - 1, ((j / 3) % 3) as i8 - 1, (j % 3) as i8 - 1);
    i += 1;
  }
  c
};

pub const UNIT_VECTORS: [Vec3; 6] =
  [Vec3(-1, 0, 0), Vec3(0, -1, 0), Vec3(0, 0, -1), Vec3(0, 0, 1), Vec3(0, 1, 0), Vec3(1, 0, 0)];

pub const MOVES: [Move; 12] = {
  let (mut moves, mut i) = ([Move { normal: Vec3(0, 0, 0), dir: 0 }; 12], 0);
  while i < 12 {
    moves[i] = Move { normal: UNIT_VECTORS[i / 2], dir: if i % 2 == 0 { -1 } else { 1 } };
    i += 1;
  }
  moves
};

pub const SOLVED_CUBE: Cube = {
  let (mut c, mut i) = ([(Vec3(0, 0, 0), Mat3::ID); 26], 0);
  while i < 26 {
    c[i] = (CUBELETS[i], Mat3::ID);
    i += 1;
  }
  c
};

// Static layer partitioning indices to avoid scanning all 26 cubelets in heuristic evaluations.
const TOP_EDGES: [usize; 4] = [5, 11, 16, 22];
const TOP_CUBELETS: [usize; 9] = [2, 5, 8, 11, 13, 16, 19, 22, 25];
const MID_CUBELETS: [usize; 17] = [1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 16, 18, 19, 21, 22, 24, 25];
const MID_LAYER: [usize; 8] = [1, 4, 7, 10, 15, 18, 21, 24];
const BOT_EDGES: [usize; 4] = [3, 9, 14, 20];
const BOT_CORNERS: [usize; 4] = [0, 6, 17, 23];

pub fn is_cubelet_solved(c: Vec3, r: Mat3) -> bool {
  (r * Mat3::diag(c)) == Mat3::diag(c)
}

pub fn is_cubelet_pos_solved(c: Vec3, r: Mat3) -> bool {
  (r * c) == c
}

const SQRT: [f32; 4] = [0.0, 1.0, std::f32::consts::SQRT_2, 1.732_050_8];

// -- Cayley Group Precomputation: L1-resident action & distance tables --------

/// Precomputed tables for the finite rotation group SO(3, ℤ) of order 24.
/// The Cayley distance graph under 12 quarter-turn moves has diameter ≤ 3.
pub struct Tables {
  pub rotations: [Mat3; 24],
  pub transition: [[[u8; 24]; 26]; 12],
  pub sqrt_dist_solved: [[f32; 24]; 26],
  pub sqrt_dist_pos: [[f32; 24]; 26],
  pub is_solved: [[bool; 24]; 26],
  pub is_pos_solved: [[bool; 24]; 26],
  pub prune_move: [[bool; 12]; 13],
}

impl Tables {
  pub fn get() -> &'static Tables {
    static INSTANCE: OnceLock<Tables> = OnceLock::new();
    INSTANCE.get_or_init(Self::init)
  }

  fn init() -> Self {
    let (mut rotations, mut i) = (vec![Mat3::ID], 0);
    while i < rotations.len() {
      for m in &MOVES {
        let nxt = m.rot_mat() * rotations[i];
        if rotations.len() < 24 && !rotations.contains(&nxt) {
          rotations.push(nxt);
        }
      }
      i += 1;
    }
    let rotations: [Mat3; 24] = rotations.try_into().unwrap();

    let (mut transition, mut rot_trans) = ([[[0u8; 24]; 26]; 12], [[0u8; 24]; 12]);
    for (m, &mv) in MOVES.iter().enumerate() {
      let (n, r_m) = (mv.normal, mv.rot_mat());
      for (r_idx, &r) in rotations.iter().enumerate() {
        let next_r = rotations.iter().position(|&x| x == r_m * r).unwrap() as u8;
        rot_trans[m][r_idx] = next_r;
        for (c_idx, &c) in CUBELETS.iter().enumerate() {
          transition[m][c_idx][r_idx] = if n.dot(r * c) > 0 { next_r } else { r_idx as u8 };
        }
      }
    }

    let compute_tables = |goal: fn(Vec3, Mat3) -> bool| {
      let mut dist = [[255u8; 24]; 26];
      for (c, &pos) in CUBELETS.iter().enumerate() {
        let (mut q, mut qh) = (Vec::with_capacity(24), 0);
        for (r, &rot) in rotations.iter().enumerate() {
          if goal(pos, rot) {
            dist[c][r] = 0;
            q.push(r);
          }
        }
        while qh < q.len() {
          let curr = q[qh];
          qh += 1;
          for &n in &rot_trans {
            let nxt = n[curr] as usize;
            if dist[c][nxt] == 255 {
              dist[c][nxt] = dist[c][curr] + 1;
              q.push(nxt);
            }
          }
        }
      }
      (
        std::array::from_fn(|c| std::array::from_fn(|r| SQRT[dist[c][r] as usize])),
        std::array::from_fn(|c| std::array::from_fn(|r| dist[c][r] == 0)),
      )
    };

    let (sqrt_dist_solved, is_solved) = compute_tables(is_cubelet_solved);
    let (sqrt_dist_pos, is_pos_solved) = compute_tables(is_cubelet_pos_solved);

    let mut prune_move = [[false; 12]; 13];
    for last in 1..=12 {
      for next in 0..12 {
        let (p, m) = (MOVES[last - 1], MOVES[next]);
        prune_move[last][next] =
          (p.normal == m.normal && p.dir == -m.dir) || (p.normal.dot(m.normal) == -1 && p.normal > m.normal);
      }
    }

    Self { rotations, transition, sqrt_dist_solved, sqrt_dist_pos, is_solved, is_pos_solved, prune_move }
  }
}

// -- Flat 26-Byte State Representation & Heuristics ---------------------------

/// Flat 26-byte cube representation. Index `i` is home position `CUBELETS[i]`,
/// and byte `self.0[i]` is the orientation index in `Tables::rotations` (0..24).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
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
    Self(std::array::from_fn(|i| {
      let &(_, r) = cube.iter().find(|(c, _)| *c == CUBELETS[i]).unwrap();
      t.rotations.iter().position(|&x| x == r).unwrap() as u8
    }))
  }

  pub fn to_cube(&self) -> Cube {
    std::array::from_fn(|i| (CUBELETS[i], Tables::get().rotations[self.0[i] as usize]))
  }

  #[inline(always)]
  pub fn apply_move(&self, m: usize) -> Self {
    Self(std::array::from_fn(|i| Tables::get().transition[m][i][self.0[i] as usize]))
  }

  #[inline(always)]
  pub fn step(&mut self, m: usize) -> usize {
    *self = self.apply_move(m);
    m
  }

  #[inline(always)]
  pub fn is_solved(&self) -> bool {
    (0..26).all(|i| Tables::get().is_solved[i][self.0[i] as usize])
  }

  #[inline(always)]
  pub fn sum_dist(&self, idx: &[usize], dist: &[[f32; 24]; 26]) -> f32 {
    idx.iter().map(|&i| dist[i][self.0[i] as usize]).sum()
  }

  #[inline(always)]
  pub fn layer_h(&self, idx: &[usize], dist: &[[f32; 24]; 26], div: f32) -> f64 {
    let s = self.sum_dist(idx, dist);
    ((s * s) / div) as f64
  }
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

// -- Pruned A* Search with Stochastic Restarts --------------------------------

pub fn astar<G: Fn(&FastCube) -> bool, H: Fn(&FastCube) -> f64>(
  start: FastCube,
  is_goal: G,
  heuristic: H,
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
          Entry::Occupied(mut e) if e.get().0 > next_cost => {
            e.insert((next_cost, mv_code));
          }
          Entry::Vacant(e) => {
            e.insert((next_cost, mv_code));
          }
          _ => continue,
        }

        if is_goal(&dst) {
          let (mut path, mut curr) = (Vec::new(), dst);
          while let Some(&(_, code @ 1..)) = visited.get(&curr) {
            path.push((code - 1) as usize);
            curr = curr.apply_move((code - 1) as usize ^ 1);
          }
          path.reverse();
          return Some((dst, path));
        }

        let hw = if random_weight > 0.0 { rng.random_gauss(1.0, random_weight).max(0.01) } else { 1.0 };
        // Non-negative IEEE-754 floats map monotonically to u64 integer bit patterns.
        frontier.push((Reverse(((next_cost as f64) + hw * heuristic(&dst)).to_bits()), next_cost, mv_code, dst));
      }
    }

    if simulated < budget || random_weight <= 0.0 {
      return None;
    }
    log(&format!("search budget of {} moves exceeded; restarting", budget));
    budget = (budget + 1).max((budget as f64 * 1.5) as usize);
  }
}

// -- Canonical Reduction Phases & 3D Solver -----------------------------------

fn top_done(t: &Tables, i: usize, c: &FastCube) -> bool {
  let cnt = |s: &[usize]| s.iter().filter(|&&j| t.is_solved[j][c.0[j] as usize]).count();
  cnt(&TOP_EDGES) >= 4.min(i + 1) && cnt(&TOP_CUBELETS) >= 9.min(i + 1) && cnt(&MID_CUBELETS) >= 17.min(i + 1)
}

fn cross_done(t: &Tables, i: usize, c: &FastCube) -> bool {
  let cnt = |s: &[usize]| s.iter().filter(|&&j| t.is_solved[j][c.0[j] as usize]).count();
  let d = Vec3(0, 0, -1);
  let ok = BOT_EDGES.iter().filter(|&&j| t.rotations[c.0[j] as usize] * d == d).count();
  top_done(t, 16, c) && ok >= 4.min(i + 1) && cnt(&BOT_EDGES) >= 4.min(i.saturating_sub(3))
}

fn corners_done(t: &Tables, i: usize, c: &FastCube) -> bool {
  let pos = BOT_CORNERS.iter().filter(|&&j| t.is_pos_solved[j][c.0[j] as usize]).count();
  cross_done(t, 7, c) && pos >= 4.min(i + 1)
}

/// Generic, monomorphized phase runner that logs progress and executes A* with random restarts.
fn solve_phase<G: Fn(usize, &FastCube) -> bool, H: Fn(usize, &FastCube) -> f64>(
  name: &str,
  total: usize,
  rw: f64,
  goal: G,
  heuristic: H,
  c: &mut FastCube,
  moves: &mut Vec<usize>,
) {
  for i in 0..total {
    log(&format!("{name} #{}", i + 1));
    let (next_c, p) = astar(*c, |x| goal(i, x), |x| heuristic(i, x), rw, 30_000).expect("phase failed");
    log(&format!("-> found solution with {} moves", p.len()));
    *c = next_c;
    moves.extend(p);
  }
  log("--------------------------------------------------");
}

/// Tests if corner `c` with rotation `r` can be oriented via quarter-turns of the bottom face.
fn corner_oriented(t: &Tables, c: usize, mut r: usize, bot: usize) -> bool {
  (0..4).any(|_| {
    let ok = t.is_solved[c][r];
    r = t.transition[bot][c][r] as usize;
    ok
  })
}

pub fn solve_endgame(cube: &mut FastCube, moves: &mut Vec<usize>) {
  let t = Tables::get();
  let find_m = |n| MOVES.iter().position(|&m| m.normal == n && m.dir == 1).unwrap();
  let (left, top, bot) = (find_m(Vec3(0, -1, 0)), find_m(Vec3(0, 0, 1)), find_m(Vec3(0, 0, -1)));
  // Sexy-move commutator [L', U', L, U, L', U', L, U] orients bottom-layer corners.
  let routine = [left ^ 1, top ^ 1, left, top, left ^ 1, top ^ 1, left, top];

  for _ in 0..4 {
    while !{
      let c = (0..26).find(|&i| t.rotations[cube.0[i] as usize] * CUBELETS[i] == Vec3(1, -1, -1)).unwrap();
      corner_oriented(t, c, cube.0[c] as usize, bot)
    } {
      moves.extend(routine.map(|m| cube.step(m)));
    }
    moves.push(cube.step(bot));
  }

  while !cube.is_solved() {
    moves.push(cube.step(bot));
  }
}

pub fn shuffle(iters: usize, seed: u64, cube: Cube) -> Cube {
  let mut rng = Rng::new(seed);
  (0..iters)
    .fold(FastCube::from_cube(&cube), |c, _| c.apply_move((rng.next_f64() * 12.0) as usize % 12))
    .to_cube()
}

pub fn solve(cube: Cube) -> Vec<Move> {
  let (t0, start_sim, mut moves) = (Instant::now(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed), Vec::new());
  let mut c = FastCube::from_cube(&cube);
  let t = Tables::get();

  let h1 = |i, c: &FastCube| {
    if i < 4 {
      c.layer_h(&TOP_EDGES, &t.sqrt_dist_solved, 4.0)
    } else if i < 9 {
      c.layer_h(&TOP_CUBELETS, &t.sqrt_dist_solved, 8.0)
    } else {
      let t_s = c.sum_dist(&TOP_CUBELETS, &t.sqrt_dist_solved);
      let m_s = c.sum_dist(&MID_LAYER, &t.sqrt_dist_solved);
      let rem = (17usize.saturating_sub(i)).max(1) as f32;
      ((t_s * t_s) / 4.0 + (m_s * m_s) / rem.max(2.0)) as f64
    }
  };
  let h_cross = |i: usize, c: &FastCube| {
    let t_s = c.sum_dist(&TOP_CUBELETS, &t.sqrt_dist_solved);
    let m_s = c.sum_dist(&MID_LAYER, &t.sqrt_dist_solved);
    if i < 4 {
      let unoriented = BOT_EDGES.iter().filter(|&&j| t.rotations[c.0[j] as usize].2 .2 != 1).count() as f32;
      let rem = (4usize.saturating_sub(i)).max(1) as f32;
      ((t_s * t_s) / 5.0 + (m_s * m_s) / 4.0 + (unoriented * unoriented) / rem.max(2.0)) as f64
    } else {
      let b_s = c.sum_dist(&BOT_EDGES, &t.sqrt_dist_solved);
      let rem = (8usize.saturating_sub(i)).max(1) as f32;
      ((t_s * t_s) / 5.0 + (m_s * m_s) / 4.0 + (b_s * b_s) / rem.max(2.0)) as f64
    }
  };
  let h_corners = |i: usize, c: &FastCube| {
    let t_s = c.sum_dist(&TOP_CUBELETS, &t.sqrt_dist_solved);
    let m_s = c.sum_dist(&MID_LAYER, &t.sqrt_dist_solved);
    let e_s = c.sum_dist(&BOT_EDGES, &t.sqrt_dist_solved);
    let c_s = c.sum_dist(&BOT_CORNERS, &t.sqrt_dist_pos);
    let rem = (4usize.saturating_sub(i)).max(1) as f32;
    ((t_s * t_s) / 6.0 + (m_s * m_s) / 6.0 + (e_s * e_s) / 6.0 + (c_s * c_s) / rem.max(2.0)) as f64
  };

  solve_phase("solving cubelet", 17, 0.25, |i, c| top_done(t, i, c), h1, &mut c, &mut moves);
  solve_phase("solving bottom cross", 8, 0.25, |i, c| cross_done(t, i, c), h_cross, &mut c, &mut moves);
  solve_phase("positioning bottom corners", 4, 0.30, |i, c| corners_done(t, i, c), h_corners, &mut c, &mut moves);
  solve_endgame(&mut c, &mut moves);

  let (elapsed, sim) = (t0.elapsed().as_secs_f64(), TOTAL_MOVES_SIMULATED.load(Ordering::Relaxed) - start_sim);
  log(&format!("Solved cube in {} moves.", moves.len()));
  log(&format!("is_cube_solved: {}", c.is_solved()));
  log(&format!("- time elapsed: {:.3} sec", elapsed));
  log(&format!("- moves simulated: {} ({:.0} moves/sec)", sim, sim as f64 / elapsed.max(0.0001)));
  moves.into_iter().map(|i| MOVES[i]).collect()
}
