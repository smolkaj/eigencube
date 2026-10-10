// eigencube.rs - Minimalistic Rubik's Cube Solver in Rust
// A Functional Pearl: Discrete 3D Euclidean space, linear algebra,
// and multi-phase A* search with move-budgeted restarts.

use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::ops::Mul;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
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
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub struct Mat3(pub Vec3, pub Vec3, pub Vec3);

impl Mat3 {
    pub const ID: Mat3 = Mat3(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));

    pub const fn diag(v: Vec3) -> Mat3 {
        Mat3(Vec3(v.0, 0, 0), Vec3(0, v.1, 0), Vec3(0, 0, v.2))
    }

    pub const fn transpose(self) -> Mat3 {
        Mat3(
            Vec3((self.0).0, (self.1).0, (self.2).0),
            Vec3((self.0).1, (self.1).1, (self.2).1),
            Vec3((self.0).2, (self.1).2, (self.2).2),
        )
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
        let c = m.transpose();
        Mat3(
            Vec3(self.0.dot(c.0) as i8, self.0.dot(c.1) as i8, self.0.dot(c.2) as i8),
            Vec3(self.1.dot(c.0) as i8, self.1.dot(c.1) as i8, self.1.dot(c.2) as i8),
            Vec3(self.2.dot(c.0) as i8, self.2.dot(c.1) as i8, self.2.dot(c.2) as i8),
        )
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Move {
    pub normal: Vec3,
    pub dir: i8,
}

impl Move {
    pub const fn invert(self) -> Move {
        Move { normal: self.normal, dir: -self.dir }
    }

    pub const fn rot_mat(self) -> Mat3 {
        let (x, y, _) = (self.normal.0, self.normal.1, self.normal.2);
        let d = self.dir;
        if x != 0 {
            Mat3(Vec3(1, 0, 0), Vec3(0, 0, d), Vec3(0, -d, 0))
        } else if y != 0 {
            Mat3(Vec3(0, 0, d), Vec3(0, 1, 0), Vec3(-d, 0, 0))
        } else {
            Mat3(Vec3(0, d, 0), Vec3(-d, 0, 0), Vec3(0, 0, 1))
        }
    }
}

pub type Cubelet = (Vec3, Mat3);
pub type Cube = [Cubelet; 26];

pub const CUBELETS: [Vec3; 26] = {
    let mut arr = [Vec3(0, 0, 0); 26];
    let coords: [i8; 3] = [-1, 0, 1];
    let mut i = 0;
    let mut xi = 0;
    while xi < 3 {
        let mut yi = 0;
        while yi < 3 {
            let mut zi = 0;
            while zi < 3 {
                let v = Vec3(coords[xi], coords[yi], coords[zi]);
                if v.norm1() > 0 {
                    arr[i] = v;
                    i += 1;
                }
                zi += 1;
            }
            yi += 1;
        }
        xi += 1;
    }
    arr
};

pub const UNIT_VECTORS: [Vec3; 6] =
    [Vec3(-1, 0, 0), Vec3(0, -1, 0), Vec3(0, 0, -1), Vec3(0, 0, 1), Vec3(0, 1, 0), Vec3(1, 0, 0)];

pub const MOVES: [Move; 12] = {
    let mut arr = [Move { normal: Vec3(0, 0, 0), dir: 0 }; 12];
    let mut i = 0;
    let mut vi = 0;
    while vi < 6 {
        arr[i] = Move { normal: UNIT_VECTORS[vi], dir: -1 };
        arr[i + 1] = Move { normal: UNIT_VECTORS[vi], dir: 1 };
        i += 2;
        vi += 1;
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
    let colors = Mat3::diag(c);
    (r * colors) == colors
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

pub fn should_prune(last_move: Option<Move>, m: Move) -> bool {
    match last_move {
        None => false,
        Some(prev) => {
            (prev.normal == m.normal && prev.dir == -m.dir)
                || (prev.normal.dot(m.normal) == -1 && prev.normal > m.normal)
        }
    }
}

pub fn log(msg: &str) {
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
    let sec = now % 60;
    let min = (now / 60) % 60;
    let hour = (now / 3600) % 24;
    println!("[{:02}:{:02}:{:02}] {}", hour, min, sec, msg);
}

// Minimal, zero-dependency Xorshift64 PRNG
pub struct Rng(pub u64);
impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
    pub fn random_gauss(&mut self, mean: f64, stdev: f64) -> f64 {
        let mut sum = 0.0;
        for _ in 0..12 {
            sum += self.next_f64();
        }
        mean + stdev * (sum - 6.0)
    }
}

// Multi-phase A* search with move-budgeted restarts (1.5x expansion)
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct VisitedNode<S> {
    pub cost: u32,
    pub prev: Option<(S, Move)>,
}

#[derive(Clone, PartialEq)]
pub struct FrontierNode<S> {
    pub prio: f64,
    pub cost: u32,
    pub last_move: Option<Move>,
    pub state: S,
}

impl<S: PartialEq> Eq for FrontierNode<S> {}
impl<S: PartialEq> Ord for FrontierNode<S> {
    fn cmp(&self, other: &Self) -> Ordering {
        other.prio.partial_cmp(&self.prio).unwrap_or(Ordering::Equal)
    }
}
impl<S: PartialEq> PartialOrd for FrontierNode<S> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn astar<S: Copy + Eq + std::hash::Hash>(
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
    let mut visited: FxHashMap<S, VisitedNode<S>> = FxHashMap::default();
    visited.reserve(65536);
    let mut rng = Rng::new(42);
    let mut budget = max_moves;

    loop {
        visited.clear();
        visited.insert(start, VisitedNode { cost: 0, prev: None });
        let mut frontier = BinaryHeap::new();
        frontier.push(FrontierNode { prio: 0.0, cost: 0, last_move: None, state: start });
        let mut simulated = 0;
        let mut goal_node = None;

        while let Some(FrontierNode { cost, last_move, state, .. }) = frontier.pop() {
            for &m in &MOVES {
                if simulated >= budget {
                    break;
                }
                if should_prune(last_move, m) {
                    continue;
                }
                TOTAL_MOVES_SIMULATED.fetch_add(1, AtomicOrdering::Relaxed);
                simulated += 1;
                let next_cost = cost + 1;
                let dst = apply_mv(m, &state);

                match visited.get(&dst) {
                    Some(v) if v.cost <= next_cost => continue,
                    _ => {
                        visited
                            .insert(dst, VisitedNode { cost: next_cost, prev: Some((state, m)) });
                        if is_goal(&dst) {
                            goal_node = Some(dst);
                            break;
                        }
                        let hw = if random_weight > 0.0 {
                            rng.random_gauss(1.0, random_weight).max(0.01)
                        } else {
                            1.0
                        };
                        let prio = (next_cost as f64) + hw * heuristic(&dst);
                        frontier.push(FrontierNode {
                            prio,
                            cost: next_cost,
                            last_move: Some(m),
                            state: dst,
                        });
                    }
                }
            }
            if goal_node.is_some() || simulated >= budget {
                break;
            }
        }

        if let Some(target) = goal_node {
            let mut path = Vec::new();
            let mut curr = target;
            while let Some(v) = visited.get(&curr) {
                if let Some((prev_state, m)) = v.prev {
                    path.push(m);
                    curr = prev_state;
                } else {
                    break;
                }
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

// Single-cubelet distance heuristics via unified A* search
type HeuristicMap = FxHashMap<(Vec3, Mat3), usize>;
static HEURISTIC_CACHE: OnceLock<(HeuristicMap, HeuristicMap)> = OnceLock::new();

fn get_heuristic_caches() -> &'static (HeuristicMap, HeuristicMap) {
    HEURISTIC_CACHE.get_or_init(|| {
        let mut solved_cache = FxHashMap::default();
        let mut pos_cache = FxHashMap::default();
        for &c in &CUBELETS {
            // Find all reachable orientations
            let mut visited = FxHashMap::default();
            let mut queue = std::collections::VecDeque::new();
            visited.insert(Mat3::ID, 0usize);
            queue.push_back(Mat3::ID);
            while let Some(r) = queue.pop_front() {
                let dist = visited[&r];
                for &m in &MOVES {
                    let next_r = m.rot_mat() * r;
                    if let std::collections::hash_map::Entry::Vacant(e) = visited.entry(next_r) {
                        e.insert(dist + 1);
                        queue.push_back(next_r);
                    }
                }
            }
            // For each reachable orientation, compute distance to solved and pos_solved via A*
            for &r in visited.keys() {
                let solved_d = astar(
                    r,
                    |&r_curr| is_cubelet_solved(c, r_curr),
                    |m, &r_curr| m.rot_mat() * r_curr,
                    |_| 0.0,
                    0.0,
                    100_000,
                )
                .map(|(_, p)| p.len())
                .unwrap_or(0);
                let pos_d = astar(
                    r,
                    |&r_curr| is_cubelet_pos_solved(c, r_curr),
                    |m, &r_curr| m.rot_mat() * r_curr,
                    |_| 0.0,
                    0.0,
                    100_000,
                )
                .map(|(_, p)| p.len())
                .unwrap_or(0);
                solved_cache.insert((c, r), solved_d);
                pos_cache.insert((c, r), pos_d);
            }
        }
        (solved_cache, pos_cache)
    })
}

pub fn min_moves_to_solved(c: Vec3, r: Mat3) -> usize {
    get_heuristic_caches().0.get(&(c, r)).copied().unwrap_or(0)
}

pub fn min_moves_to_pos(c: Vec3, r: Mat3) -> usize {
    get_heuristic_caches().1.get(&(c, r)).copied().unwrap_or(0)
}

pub fn norm_p05(cube: &Cube, cond: impl Fn(Vec3) -> bool, f: impl Fn(Vec3, Mat3) -> usize) -> f64 {
    let sum: f64 =
        cube.iter().filter(|&&(c, _)| cond(c)).map(|&(c, r)| (f(c, r) as f64).sqrt()).sum();
    sum * sum
}

pub fn top_layer_heuristic(cube: &Cube) -> f64 {
    norm_p05(cube, |c| c.2 == 1, min_moves_to_solved) / 8.0
}

pub fn middle_layer_heuristic(cube: &Cube) -> f64 {
    norm_p05(cube, |c| c.2 >= 0, min_moves_to_solved) / 4.0
}

pub fn bottom_layer_edge_heuristic(cube: &Cube) -> f64 {
    norm_p05(cube, |c| !(c.2 == -1 && c.norm1() == 3), min_moves_to_solved) / 3.0
}

pub fn bottom_layer_corner_heuristic(cube: &Cube) -> f64 {
    (norm_p05(cube, |c| c.2 == 1, min_moves_to_solved) / 5.0)
        + (norm_p05(cube, |c| c.2 == 0, min_moves_to_solved) / 3.0)
        + (norm_p05(
            cube,
            |c| c.2 == -1,
            |c, r| {
                if c.norm1() == 3 {
                    min_moves_to_pos(c, r)
                } else {
                    min_moves_to_solved(c, r)
                }
            },
        ) / 8.0)
}

pub fn count_solved(cube: &Cube, f: impl Fn(Vec3) -> bool) -> usize {
    cube.iter().filter(|&&(c, r)| f(c) && is_cubelet_solved(c, r)).count()
}

pub fn count_bottom_edges_positioned(cube: &Cube) -> usize {
    cube.iter()
        .filter(|&&(c, r)| c.2 == -1 && c.norm1() == 2 && (r * Vec3(0, 0, -1)) == Vec3(0, 0, -1))
        .count()
}

pub fn count_bottom_corners_positioned(cube: &Cube) -> usize {
    cube.iter()
        .filter(|&&(c, r)| c.2 == -1 && c.norm1() == 3 && is_cubelet_pos_solved(c, r))
        .count()
}

pub fn is_top_edge(v: Vec3) -> bool {
    v.2 == 1 && v.norm1() == 2
}
pub fn is_bottom_edge(v: Vec3) -> bool {
    v.2 == -1 && v.norm1() == 2
}

fn solve_layer(
    name: &str,
    total: usize,
    is_goal: impl Fn(usize, &Cube) -> bool,
    heuristic: impl Fn(usize) -> fn(&Cube) -> f64,
    random_weight: f64,
    mut cube: Cube,
) -> (Cube, Vec<Move>) {
    let mut all_moves = Vec::new();
    for i in 0..total {
        log(&format!("{} #{}", name, i + 1));
        let h = heuristic(i);
        match astar(cube, |c| is_goal(i, c), apply_move, h, random_weight, 100_000) {
            Some((next_c, mvs)) => {
                log(&format!("-> found solution with {} moves", mvs.len()));
                cube = next_c;
                all_moves.extend(mvs);
            }
            None => panic!("Failed {}", name),
        }
    }
    (cube, all_moves)
}

pub fn bottom_left_front_corner(cube: &Cube) -> Cubelet {
    *cube.iter().find(|&&(c, r)| (r * c) == Vec3(1, -1, -1)).expect("corner not found")
}

// Endgame: orient bottom corners using (R' D' R D) * 2/4 and align bottom face
pub fn solve_endgame(cube: Cube) -> (Cube, Vec<Move>) {
    let left = Move { normal: Vec3(0, -1, 0), dir: 1 };
    let top = Move { normal: Vec3(0, 0, 1), dir: 1 };
    let bottom = Move { normal: Vec3(0, 0, -1), dir: 1 };
    let cycle = [left.invert(), top.invert(), left, top];
    let routine = [cycle[0], cycle[1], cycle[2], cycle[3], cycle[0], cycle[1], cycle[2], cycle[3]];

    let mut state_cube = cube;
    let mut solution = Vec::new();

    let apply = |m: Move, c: &mut Cube, sol: &mut Vec<Move>| {
        *c = apply_move(m, c);
        sol.push(m);
    };

    let is_corner_oriented = |c: &Cube| -> bool {
        let (c_orig, mut r) = bottom_left_front_corner(c);
        let r_prime = bottom.rot_mat();
        for _ in 0..4 {
            if is_cubelet_solved(c_orig, r) {
                return true;
            }
            r = r_prime * r;
        }
        false
    };

    for _ in 0..4 {
        while !is_corner_oriented(&state_cube) {
            for &m in &routine {
                apply(m, &mut state_cube, &mut solution);
            }
        }
        apply(bottom, &mut state_cube, &mut solution);
    }

    while !is_cube_solved(&state_cube) {
        apply(bottom, &mut state_cube, &mut solution);
    }

    (state_cube, solution)
}

pub fn shuffle(iters: usize, seed: u64, mut cube: Cube) -> Cube {
    let mut rng = Rng::new(seed);
    for _ in 0..iters {
        let m = MOVES[rng.next_u64() as usize % MOVES.len()];
        cube = apply_move(m, &cube);
    }
    cube
}

// Full 3-phase human solver: top layer -> middle edges -> bottom layer & endgame
pub fn solve(cube: Cube) -> Vec<Move> {
    let t0 = Instant::now();
    let start_sim = TOTAL_MOVES_SIMULATED.load(AtomicOrdering::Relaxed);

    let (c1, s1) = solve_layer(
        "solving cubelet",
        17,
        |i, c| {
            count_solved(c, is_top_edge) >= 4.min(i + 1)
                && count_solved(c, |v| v.2 == 1) >= 9.min(i + 1)
                && count_solved(c, |v| v.2 >= 0) >= 17.min(i + 1)
        },
        |i| {
            if i < 9 {
                top_layer_heuristic
            } else {
                middle_layer_heuristic
            }
        },
        0.25,
        cube,
    );
    log("--------------------------------------------------");

    let (c2, s2) = solve_layer(
        "solving bottom cross",
        8,
        |i, c| {
            count_solved(c, |v| v.2 >= 0) == 17
                && count_bottom_edges_positioned(c) >= 4.min(i + 1)
                && count_solved(c, is_bottom_edge) >= 4.min(i.saturating_sub(3))
        },
        |_| bottom_layer_edge_heuristic,
        0.25,
        c1,
    );
    log("--------------------------------------------------");

    let (c3, s3) = solve_layer(
        "positioning bottom corners",
        4,
        |i, c| {
            count_solved(c, |v| v.2 >= 0) == 17
                && count_solved(c, is_bottom_edge) == 4
                && count_bottom_corners_positioned(c) >= 4.min(i + 1)
        },
        |_| bottom_layer_corner_heuristic,
        0.30,
        c2,
    );
    log("--------------------------------------------------");

    let (c4, s4) = solve_endgame(c3);
    let mut moves = s1;
    moves.extend(s2);
    moves.extend(s3);
    moves.extend(s4);

    let elapsed = t0.elapsed().as_secs_f64();
    let moves_simulated = TOTAL_MOVES_SIMULATED.load(AtomicOrdering::Relaxed) - start_sim;
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
