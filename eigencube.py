""" A minimalstic Rubik's cube solver.

We model a Rubik's cube as a discrete object in 3-dimensional Euclidean space,
centered at the origin (0, 0, 0).
Each *cubelet* has coordinates (x, y, z) in {-1, 0, 1}^3.
Each *move* is a 90 degree hyperplane rotation, with the hyperplane given
by a standard unit vector or its opposite.
"""

import numpy as np
import heapq
import random
import functools
import itertools
from datetime import datetime
from dataclasses import dataclass, field
from typing import Any

# Whether or not to use randomization during search.
RANDOMIZE_SEARCH = True

# Remember start up time for stats.
STARTUP_TIME = datetime.now()

# Redefine print to include timestamps.
_print = print
def print(*args, **kw):
  _print("[%s]" % (datetime.now().strftime('%H:%M:%S')), *args, **kw)

# Returns the 1-norm of a vector.
def norm1(v): return abs(v[0]) + abs(v[1]) + abs(v[2])

# Returns a 2-dimensional matrix as a tuple.
def tupled(np_mat): return tuple(tuple(int(x) for x in row) for row in np_mat)

# We are interested in x, y, z coordinates ranging over {-1, 0, 1}.
crange = [-1, 0, 1]
vectors = tuple((x,y,z) for x in crange for y in crange for z in crange)
unit_vectors = [v for v in vectors if norm1(v) == 1]

# A cube is encoded as a finite map from cubelets to 3-by-3 rotation matrices
# `r`, where each cubelet is encoded as the 3D vector `c` that indicates the
# position of the cubelet in the solved cube.
# The cubelets current position is given by the matrix-vector product `c*r`.
solved_cube = tuple((c, tupled(np.identity(3))) for c in vectors if any(c))
NUM_CUBELETS = len(solved_cube)

# A move is a clockwise or counterclockwise 90 degree rotation of the
# slice pointed at by a unit vector.
moves = [(v, direction) for v in unit_vectors for direction in [-1, 1]]

# Each color is encoded by the unit vector corresponding to the direction that
# faces of that color point to in a solved cube.
color_names = {
  (+1, 0, 0): "GREEN",   # front
  (0, +1, 0): "RED",     # right
  (0, 0, +1): "WHITE",   # top
  (-1, 0, 0): "BLUE",    # back
  (0, -1, 0): "ORANGE",  # left
  (0, 0, -1): "YELLOW",  # bottom
}
assert all(v in color_names for v in unit_vectors)

cubelet_types = ("hidden", "center", "edge", "corner")
def describe_cubelet_type(cubelet): return cubelet_types[norm1(cubelet)]

def describe_vector(vector):
  x, y, z = vector
  descriptions = []
  descriptions += ["top"] if z == 1 else ["bottom"] if z == -1 else []
  descriptions += ["right"] if y == 1 else ["left"] if y == -1 else []
  descriptions += ["front"] if x == 1 else ["back"] if x == -1 else []
  assert descriptions
  return "-".join(descriptions)

def describe_config(cubelet, rotation):
  colors = np.diag(cubelet)
  sticker_directions = rotation @ colors
  return ", ".join(sorted(
     describe_vector(direction) + ": " + color_names[tuple(color)]
     for color, direction in zip(colors.T, sticker_directions.T)
     if any(color)
  ))

@functools.cache
def position(cubelet, rotation): return tuple(np.matmul(rotation, cubelet))

def describe_cubelet(cubelet, rotation):
  return "%s %s: %s" % (
    describe_cubelet_type(cubelet),
    describe_vector(position(cubelet, rotation)),
    describe_config(cubelet, rotation),
  )

def describe_cube(cube):
  return "\n".join(sorted(describe_cubelet(*c) for c in cube))

def describe_move(move):
  v, direction = move
  return "%s rotation of %s slice" % (
      "clockwise" if direction == 1 else "counterclockwise",
      describe_vector(v),
  )

def inverse_move(move):
  v, direction = move
  return (v, -direction)

@functools.cache
def rotation_matrix(move):
  v, direction = move
  # The rotational axis is the dimension into which `v` is pointing.
  fixed_dim = next(i for i in range(3) if v[i])
  # The rotation takes place in the other two dimensions.
  r1, r2 = (i for i in range(3) if i != fixed_dim)

  M = np.zeros([3, 3])
  M[fixed_dim, fixed_dim] = 1
  # The 90 degree rotation matrix [[0, 1], [-1,0]], adjusted by direction.
  # https://en.wikipedia.org/wiki/Rotation_matrix#Common_rotations
  M[r1, r2], M[r2, r1] = direction, -direction
  return M

@functools.cache
def apply_move_to_cubelet_rotation(move, cubelet, rotation):
  v, direction = move
  move_applies = np.dot(v, position(cubelet, rotation)) > 0
  return tupled(rotation_matrix(move) @ rotation) if move_applies else rotation

def apply_move_to_cube(move, cube):
  return tuple(
    (cubelet, apply_move_to_cubelet_rotation(move, cubelet, rotation))
    for cubelet, rotation in cube
  )

def shuffle(cube, iterations=100_000, seed=42):
  if seed is not None: random.seed(seed)
  for _ in range(iterations):
    move = moves[random.randrange(len(moves))]
    cube = apply_move_to_cube(move, cube)
  return cube

def run_tests():
  # Check that every move is a permutation with cycle length 4.
  for move in moves:
    cubes = [solved_cube]
    for _ in range(3): cubes.append(apply_move_to_cube(move, cubes[-1]))
    assert len(set(cubes)) == 4
    assert apply_move_to_cube(move, cubes[-1]) == cubes[0]


@dataclass(order=True, frozen=True)
class PrioritizedItem:
  item: Any = field(compare=False)
  priority: int

# A search step is a sequence of moves: a single move, or a learned macro.
single_move_steps = tuple((move,) for move in moves)

def inverse_step(step): return tuple(inverse_move(move) for move in reversed(step))

def astar(start, is_goal, apply_step, heuristic = lambda _: 0,
          get_steps = lambda _: single_move_steps, random_weight=0, max_expansions=10_000):
  if is_goal(start): return (start, ())
  budget = max_expansions

  def reconstruct_solution(dst):
    path, current = [], dst
    while current in came_from:
      src, step = came_from[current]
      path.append(step)
      current = src
    return (dst, tuple(reversed(path)))

  while True:
    frontier = [PrioritizedItem(start, 0)]
    came_from, cost_so_far = {}, { start : 0 }
    expansions = 0

    # The budget counts expanded states, not generated ones, so that it measures
    # search progress independently of how many steps are known.
    while frontier and expansions < budget:
      expansions += 1
      src = heapq.heappop(frontier).item
      last_step = came_from[src][1] if src in came_from else None
      for step in get_steps(src):
        if last_step:
          # Never immediately undo the step just taken.
          if step == inverse_step(last_step): continue
          # Opposite face moves commute; prune duplicate branches by enforcing canonical order.
          if len(step) == len(last_step) == 1:
            (last_axis, _), (axis, _) = last_step[0], step[0]
            if last_axis > axis and tuple(-x for x in last_axis) == axis: continue
        dst, cost = apply_step(step, src), cost_so_far[src] + len(step)
        if dst in cost_so_far and cost_so_far[dst] <= cost: continue
        cost_so_far[dst], came_from[dst] = cost, (src, step)
        if is_goal(dst): return reconstruct_solution(dst)
        h_weight = random.gauss(1, random_weight) if RANDOMIZE_SEARCH else 1
        priority = cost + h_weight * heuristic(dst)
        heapq.heappush(frontier, PrioritizedItem(dst, priority))
    if not frontier or random_weight == 0 or not RANDOMIZE_SEARCH:
      return None
    print("search budget of %d expansions exceeded; restarting" % budget)
    budget = max(int(1.5 * budget), budget + 1)

@functools.cache
def is_cubelet_solved(cubelet, rotation):
  colors = np.diag(cubelet)
  sticker_directions = rotation @ colors
  return np.array_equal(sticker_directions, colors)

def is_cube_solved(cube): return all(is_cubelet_solved(c, r) for c, r in cube)

def num_solved_with_criterion(cube, criterion):
  return sum(is_cubelet_solved(c, r) for c,r in cube if criterion(c))

# Rotates an isolated cubelet by a single-move step, as if no other cubelets were in the way.
def rotate_freely(step, rotation):
  (move,) = step
  return tupled(rotation_matrix(move) @ rotation)

@functools.cache
def min_moves_to_solved(cubelet, rotation):
  def is_dst(r): return is_cubelet_solved(cubelet, r)
  _, path = astar(rotation, is_dst, rotate_freely)
  return len(path)

@functools.cache
def min_moves_to_position(cubelet, rotation):
  def is_dst(r): return position(cubelet, r) == cubelet
  _, path = astar(rotation, is_dst, rotate_freely)
  return len(path)

def top_layer_heuristic(cube):
  p, n = 0.5, 8
  d = sum(min_moves_to_solved(c, r)**p for c, r in cube if c[2] == 1) ** (1/p)
  return d/n

def middle_layer_heuristic(cube):
  p, n = 0.5, 4
  d = sum(min_moves_to_solved(c, r)**p for c, r in cube if c[2] >= 0) ** (1/p)
  return d/n

def bottom_layer_edge_heuristic(cube):
  p, n = 0.5, 3
  d = sum(min_moves_to_solved(c, r)**p for c, r in cube
          if not (c[2] == -1 and norm1(c) == 3)) ** (1/p)
  return d/n

def bottom_layer_corner_heuristic(cube):
  p, n1, n2, n3 = 0.5, 5, 3, 8
  d1 = sum(min_moves_to_solved(c, r)**p for c, r in cube if c[2] == 1) ** (1/p)
  d2 = sum(min_moves_to_solved(c, r)**p for c, r in cube if c[2] == 0) ** (1/p)
  d3 = sum((min_moves_to_position(c, r) if norm1(c) == 3 else min_moves_to_solved(c, r))**p
           for c, r in cube if c[2] == -1) ** (1/p)
  return d1/n1 + d2/n2 + d3/n3

def is_top_edge(cubelet): return cubelet[2] == 1 and norm1(cubelet) == 2
def is_top_cubelet(cubelet): return cubelet[2] == 1
def is_top_or_middle_cubelet(cubelet): return cubelet[2] >= 0

# The steps the solver searches over: for each effect on the solved cube, the
# shortest known move sequence achieving it. Beyond single moves, the solver
# learns macros from its own searches: sequences that disturb only a few
# cubelets (commutators, as it turns out). Nothing about them is hard-coded.
MAX_MACRO_DISTURBANCE = 12
steps = {}

# The 48 symmetries of the cube: signed permutation matrices.
symmetries = [np.diag(signs) @ np.identity(3)[list(permutation)]
              for permutation in itertools.permutations(range(3))
              for signs in itertools.product((-1, 1), repeat=3)]

@functools.cache
def apply_step_to_cubelet_rotation(step, cubelet, rotation):
  for move in step: rotation = apply_move_to_cubelet_rotation(move, cubelet, rotation)
  return rotation

def apply_step_to_cube(step, cube):
  return tuple((c, apply_step_to_cubelet_rotation(step, c, r)) for c, r in cube)

def num_disturbed(cube): return sum(not is_cubelet_solved(c, r) for c, r in cube)

def symmetric_move(symmetry, move):
  # Conjugating a slice rotation by a symmetry yields the rotation of the image slice.
  rotation = symmetry @ rotation_matrix(move) @ symmetry.T
  return next(m for m in moves if np.array_equal(rotation_matrix(m), rotation) and
              m[0] == tuple(int(x) for x in symmetry @ move[0]))

def learn_step(sequence):
  if not 0 < num_disturbed(apply_step_to_cube(sequence, solved_cube)) <= MAX_MACRO_DISTURBANCE: return
  for symmetry, variant in itertools.product(symmetries, (sequence, inverse_step(sequence))):
    step = tuple(symmetric_move(symmetry, move) for move in variant)
    effect = apply_step_to_cube(step, solved_cube)
    if effect not in steps or len(step) < len(steps[effect]): steps[effect] = step

for step in single_move_steps: learn_step(step)

def is_bottom_edge(cubelet): return cubelet[2] == -1 and norm1(cubelet) == 2
def is_bottom_corner(cubelet): return cubelet[2] == -1 and norm1(cubelet) == 3
def has_yellow_bottom(cubelet, rotation):
  return cubelet[2] == -1 and position((0, 0, -1), rotation) == (0, 0, -1)
def is_in_right_place(c, r): return position(c, r) == c
def num_bottom_edges_positioned(cube):
  return sum(is_bottom_edge(c) and has_yellow_bottom(c, r) for c, r in cube)
def num_bottom_corners_positioned(cube):
  return sum(is_bottom_corner(c) and is_in_right_place(c, r) for c, r in cube)

# Each phase solves a few more cubelets per step, keeping the ones already solved.
def top_layer_goal(cube, i): return (
  num_solved_with_criterion(cube, is_top_edge) >= min(4, i + 1) and
  num_solved_with_criterion(cube, is_top_cubelet) >= i + 1)
def middle_layer_goal(cube, i): return (
  num_solved_with_criterion(cube, is_top_cubelet) == 9 and
  num_solved_with_criterion(cube, is_top_or_middle_cubelet) >= 10 + i)
def bottom_cross_goal(cube, i): return (
  num_solved_with_criterion(cube, is_top_or_middle_cubelet) == 17 and
  num_bottom_edges_positioned(cube) >= min(4, i + 1) and
  num_solved_with_criterion(cube, is_bottom_edge) >= min(4, i - 3))
def bottom_corner_positions_goal(cube, i): return (
  num_solved_with_criterion(cube, is_top_or_middle_cubelet) == 17 and
  num_solved_with_criterion(cube, is_bottom_edge) == 4 and
  num_bottom_corners_positioned(cube) >= i + 1)

phases = (  # (description, number of steps, goal, heuristic, random weight)
  ("solving top layer", 9, top_layer_goal, top_layer_heuristic, 0.25),
  ("solving middle layer", 8, middle_layer_goal, middle_layer_heuristic, 0.25),
  ("solving bottom cross", 8, bottom_cross_goal, bottom_layer_edge_heuristic, 0.25),
  ("positioning bottom corners", 4, bottom_corner_positions_goal, bottom_layer_corner_heuristic, 0.3),
)

# Twists each bottom corner in place with a fixed routine, then aligns the bottom layer.
def solve_endgame(cube, report_progress_callback):
  left, top, bottom = ((0, -1, 0), 1), ((0, 0, 1), 1), ((0, 0, -1), 1)
  twist = 2 * (inverse_move(left), inverse_move(top), left, top)
  solution = ()

  def is_front_corner_oriented(cube):
    c, r = next((c, r) for c, r in cube if position(c, r) == (1, -1, -1))
    return any(is_cubelet_solved(c, apply_step_to_cubelet_rotation(k * (bottom,), c, r)) for k in range(4))

  for _ in range(4):
    report_progress_callback(cube)
    while not is_front_corner_oriented(cube):
      cube, solution = apply_step_to_cube(twist, cube), solution + twist
    cube, solution = apply_move_to_cube(bottom, cube), solution + (bottom,)
  while not is_cube_solved(cube):
    cube, solution = apply_move_to_cube(bottom, cube), solution + (bottom,)
  return (cube, solution)

def solve(cube, report_progress_callback=lambda cube: None):
  solution = ()
  for description, num_steps, goal, heuristic, random_weight in phases:
    for i in range(num_steps):
      report_progress_callback(cube)
      print("%s #%d" % (description, i + 1))
      known_steps = tuple(steps.values())
      cube, path = astar(cube, lambda cube: goal(cube, i), apply_step_to_cube, heuristic,
                         lambda _: known_steps, random_weight)
      next_moves = sum(path, ())
      print("-> found solution with %d moves (%d steps known)" % (len(next_moves), len(steps)))
      learn_step(next_moves)
      solution += next_moves
    print(50 * "-")
  cube, endgame_moves = solve_endgame(cube, report_progress_callback)
  solution += endgame_moves
  print("Solved cube in %d moves." % len(solution))
  print("is_cube_solved: ", is_cube_solved(cube))
  return solution

def print_stats():
  secs_elapsed = (datetime.now() - STARTUP_TIME).total_seconds()
  cache_info = apply_step_to_cubelet_rotation.cache_info()
  steps_simulated = (cache_info.hits + cache_info.misses) / NUM_CUBELETS
  print("- time elapsed: %.1f sec" % secs_elapsed)
  print("- steps simulated: %d (%.0f steps/sec) " % (steps_simulated, steps_simulated / secs_elapsed))
  cache_info = min_moves_to_solved.cache_info()
  print("- min moves to solved calculations: ", cache_info.hits + cache_info.misses)

if __name__ == "__main__":
  import sys
  run_tests()
  arg = sys.argv[1] if len(sys.argv) > 1 else "42"
  if arg in ("-h", "--help"):
    print("Usage: python eigencube.py [seed | --benchmark]")
  elif arg == "--benchmark":
    for seed in range(100):
      print("== SEED:", seed, "==========================================")
      solve(shuffle(solved_cube, iterations=100_000, seed=seed))
      print_stats()
  else:
    seed = int(arg)
    print("Solving scrambled cube (seed=%d)..." % seed)
    solve(shuffle(solved_cube, iterations=100_000, seed=seed))
    print_stats()
