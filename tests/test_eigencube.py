import unittest
import os
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from PIL import Image

from eigencube import (
    solved_cube,
    moves,
    unit_vectors,
    color_names,
    apply_move_to_cube,
    is_cube_solved,
    solve,
    astar,
    apply_step_to_cube,
    norm1,
    rotation_matrix,
    inverse_move,
    NUM_CUBELETS,
)


class TestEigencube(unittest.TestCase):
    def test_solved_cube_invariants(self):
        """Verify structural properties and L1 norm classifications of the solved cube."""
        self.assertEqual(len(solved_cube), NUM_CUBELETS)
        self.assertEqual(NUM_CUBELETS, 26)
        self.assertTrue(is_cube_solved(solved_cube))
        self.assertEqual(len(moves), 12)
        self.assertEqual(len(unit_vectors), 6)
        self.assertTrue(all(v in color_names for v in unit_vectors))

        # Check L1 norm classifications: 6 centers, 12 edges, 8 corners
        centers = [c for c, _ in solved_cube if norm1(c) == 1]
        edges = [c for c, _ in solved_cube if norm1(c) == 2]
        corners = [c for c, _ in solved_cube if norm1(c) == 3]
        interiors = [c for c, _ in solved_cube if norm1(c) == 0]

        self.assertEqual(len(centers), 6)
        self.assertEqual(len(edges), 12)
        self.assertEqual(len(corners), 8)
        self.assertEqual(len(interiors), 0)

        # Coordinate domain invariant: {-1, 0, 1}^3
        for c, r in solved_cube:
            self.assertTrue(all(x in (-1, 0, 1) for x in c))
            self.assertEqual(r, ((1, 0, 0), (0, 1, 0), (0, 0, 1)))

    def test_rotation_matrix_properties(self):
        """Verify that all rotation matrices are orthogonal and orientation-preserving (det(M) = 1 and M^T M = I)."""
        identity = np.eye(3)
        for move in moves:
            v, _ = move
            M = rotation_matrix(move)
            # Orthogonality: M^T @ M = I
            self.assertTrue(np.allclose(M.T @ M, identity), f"Move {move} is not orthogonal")
            # Orientation preserving: det(M) = +1 (SO(3))
            self.assertTrue(np.isclose(np.linalg.det(M), 1.0), f"det(M) != 1 for move {move}")
            # Rotational axis preservation: M @ v = v
            self.assertTrue(np.allclose(M @ np.array(v), np.array(v)), f"Axis not fixed for move {move}")

    def test_single_move_order_4(self):
        """Every 90-degree face turn must have order 4 (cycle of 4 returns to identity)."""
        for move in moves:
            cube = solved_cube
            states = [cube]
            for _ in range(4):
                cube = apply_move_to_cube(move, cube)
                states.append(cube)
            # Distinct intermediate states
            self.assertEqual(len(set(states[:4])), 4, f"Move {move} did not produce 4 distinct states")
            # 4th rotation returns to initial state
            self.assertEqual(cube, solved_cube, f"Move {move} x 4 did not return to solved state")

    def test_inverse_move_cancellation(self):
        """Applying a move and its inverse should be the identity."""
        for move in moves:
            inv = inverse_move(move)
            self.assertIn(inv, moves)
            cube = apply_move_to_cube(move, solved_cube)
            self.assertNotEqual(cube, solved_cube)
            cube_restored = apply_move_to_cube(inv, cube)
            self.assertEqual(cube_restored, solved_cube)

    def test_sexy_move_order_6(self):
        """The 'sexy move' (R U R' U') repeated 6 times returns the cube to its original state."""
        # Find R (Right: +y) and U (Up/Top: +z)
        r_cw = ((0, 1, 0), 1)
        r_ccw = ((0, 1, 0), -1)
        u_cw = ((0, 0, 1), 1)
        u_ccw = ((0, 0, 1), -1)

        sexy_move = [r_cw, u_cw, r_ccw, u_ccw]

        cube = solved_cube
        for iteration in range(1, 7):
            for m in sexy_move:
                cube = apply_move_to_cube(m, cube)
            if iteration < 6:
                self.assertFalse(
                    is_cube_solved(cube),
                    f"Cube solved prematurely at iteration {iteration}",
                )

        self.assertEqual(cube, solved_cube)
        self.assertTrue(is_cube_solved(cube))

    def test_solve_simple_scramble(self):
        """Solver should successfully solve a shallowly scrambled cube."""
        # Apply 2 known moves
        scramble_moves = [((1, 0, 0), 1), ((0, 1, 0), -1)]
        cube = solved_cube
        for m in scramble_moves:
            cube = apply_move_to_cube(m, cube)
        self.assertFalse(is_cube_solved(cube))

        # Solve
        solution = solve(cube)
        for m in solution:
            cube = apply_move_to_cube(m, cube)

        self.assertTrue(is_cube_solved(cube))

    def test_shuffle_reproducibility(self):
        """Shuffle with the same seed must produce identical cube states."""
        from eigencube import shuffle
        cube_a = shuffle(solved_cube, iterations=20, seed=42)
        cube_b = shuffle(solved_cube, iterations=20, seed=42)
        cube_c = shuffle(solved_cube, iterations=20, seed=99)
        self.assertEqual(cube_a, cube_b)
        self.assertNotEqual(cube_a, cube_c)
        self.assertFalse(is_cube_solved(cube_a))

    def test_every_seed_reproduces_its_scramble(self):
        """Every seed, including falsy 0, must pin the scramble regardless of prior random state."""
        import random
        from eigencube import shuffle
        for seed in range(5):
            cube_a = shuffle(solved_cube, iterations=20, seed=seed)
            random.random()  # Advance the global generator between the two shuffles.
            cube_b = shuffle(solved_cube, iterations=20, seed=seed)
            self.assertEqual(cube_a, cube_b, f"seed={seed} is not reproducible")

    def test_descriptions(self):
        """Verify describe_vector, describe_move, and describe_cubelet_type."""
        from eigencube import describe_vector, describe_move, describe_cubelet_type
        self.assertEqual(describe_vector((1, 0, 0)), "front")
        self.assertEqual(describe_vector((0, 1, 1)), "top-right")
        self.assertEqual(describe_vector((-1, -1, -1)), "bottom-left-back")
        self.assertEqual(describe_cubelet_type((0, 0, 1)), "center")
        self.assertEqual(describe_cubelet_type((1, 1, 0)), "edge")
        self.assertEqual(describe_cubelet_type((1, 1, 1)), "corner")
        self.assertEqual(
            describe_move(((1, 0, 0), 1)),
            "clockwise rotation of front slice",
        )
        self.assertEqual(
            describe_move(((0, 0, -1), -1)),
            "counterclockwise rotation of bottom slice",
        )

    def test_astar_budget_exhaustion(self):
        """When random_weight=0, exceeding max_expansions returns None instead of searching indefinitely."""
        # Scramble with 3 moves
        cube = solved_cube
        scramble = [((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 1)]
        for m in scramble:
            cube = apply_move_to_cube(m, cube)

        # Budget of 2 expansions is insufficient to solve a 3-move scramble
        res = astar(cube, is_cube_solved, apply_step_to_cube, random_weight=0, max_expansions=2)
        self.assertIsNone(res)

        # Sufficient budget succeeds
        res = astar(cube, is_cube_solved, apply_step_to_cube, random_weight=0, max_expansions=5000)
        self.assertIsNotNone(res)
        dst, path = res
        self.assertTrue(is_cube_solved(dst))
        self.assertEqual(len(path), 3)

    def test_astar_restart_expansion(self):
        """When random_weight > 0 and budget is tight, A* restarts with 1.5x budget and finds goal."""
        cube = solved_cube
        scramble = [((1, 0, 0), 1), ((0, 1, 0), -1)]
        for m in scramble:
            cube = apply_move_to_cube(m, cube)

        # Initial budget of 5 expansions will trigger restarts but budget expands by 1.5x until solved
        res = astar(cube, is_cube_solved, apply_step_to_cube, random_weight=0.25, max_expansions=5)
        self.assertIsNotNone(res)
        dst, path = res
        self.assertTrue(is_cube_solved(dst))
        self.assertLessEqual(len(path), 2)

    def test_astar_unreachable_frontier_exhaustion(self):
        """When the goal is unreachable and the frontier empties, astar returns None immediately without restarting."""
        # 1. Immediate exhaustion on 0-transition graph
        res = astar(
            0,
            lambda x: x == 5,
            lambda m, s: s,
            get_steps=lambda s: [],
            random_weight=0.25,
            max_expansions=1000,
        )
        self.assertIsNone(res)

        # 2. Finite 3-state cyclic component {0, 1, 2} with unreachable goal 99
        # Expansions will reach max_expansions=2, but once frontier empties it must return None
        dummy_move = ((1, 0, 0), 1)
        res_cyclic = astar(
            0,
            lambda x: x == 99,
            lambda m, s: (s + 1) % 3,
            get_steps=lambda s: [(dummy_move,)],
            random_weight=0.25,
            max_expansions=2,
        )
        self.assertIsNone(res_cyclic)

    def test_headless_gui_render(self):
        """Verify that eigencube_gui renders a frame headlessly without error."""
        import eigencube_gui

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "preview.png")
            result = eigencube_gui.render_frame_to_image(solved_cube, out_path)
            self.assertEqual(result, out_path)
            self.assertTrue(os.path.exists(out_path))
            self.assertGreater(os.path.getsize(out_path), 0)

            # Test solving progress frame rendering
            solving_path = os.path.join(tmp_dir, "solving.png")
            result_solving = eigencube_gui.render_frame_to_image(solved_cube, solving_path, solving_cube=solved_cube)
            self.assertEqual(result_solving, solving_path)
            self.assertTrue(os.path.exists(solving_path))
            self.assertGreater(os.path.getsize(solving_path), 0)

    def test_background_solve_thread_lifecycle(self):
        """Verify solver runs cleanly in background daemon thread with atomic progress updates."""
        import threading
        # 2-move shallow scramble
        scramble_moves = [((1, 0, 0), 1), ((0, 1, 0), -1)]
        cube = solved_cube
        for m in scramble_moves:
            cube = apply_move_to_cube(m, cube)

        progress_history = []
        result_holder = []

        def worker(c):
            def progress(pc):
                progress_history.append(pc)
            sol = solve(c, progress)
            result_holder.append(sol)

        t = threading.Thread(target=worker, args=(cube,), daemon=True)
        t.start()
        t.join(timeout=10.0)

        self.assertFalse(t.is_alive(), "Worker thread timed out")
        self.assertEqual(len(result_holder), 1)
        self.assertGreater(len(progress_history), 0)

        # Verify produced solution
        c = cube
        for m in result_holder[0]:
            c = apply_move_to_cube(m, c)
        self.assertTrue(is_cube_solved(c))

    def test_opposite_face_moves_commute(self):
        """Opposite face moves act on disjoint slices and commute: m1 * m2 == m2 * m1."""
        for m1 in moves:
            v1, _ = m1
            opposite_v = tuple(-x for x in v1)
            for d2 in [-1, 1]:
                m2 = (opposite_v, d2)
                # Apply m1 then m2
                s1 = apply_move_to_cube(m2, apply_move_to_cube(m1, solved_cube))
                # Apply m2 then m1
                s2 = apply_move_to_cube(m1, apply_move_to_cube(m2, solved_cube))
                self.assertEqual(s1, s2, f"Opposite moves {m1} and {m2} did not commute")

    def test_min_moves_to_position(self):
        """Verify min_moves_to_position ignores cubelet orientation."""
        from eigencube import min_moves_to_position, min_moves_to_solved, position
        corner = (1, 1, -1)
        # Identity rotation: 0 moves to solved and position
        identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        self.assertEqual(min_moves_to_position(corner, identity), 0)
        self.assertEqual(min_moves_to_solved(corner, identity), 0)

        # Find a rotation where corner is in home place but twisted
        from eigencube import shuffle
        for seed in range(50):
            cube = shuffle(solved_cube, iterations=20, seed=seed)
            for c, r in cube:
                if c == corner and position(c, r) == corner and r != identity:
                    # Corner is positioned but twisted
                    self.assertEqual(min_moves_to_position(c, r), 0)
                    self.assertGreater(min_moves_to_solved(c, r), 0)
                    return

    def test_headless_gui_render_with_solution(self):
        """Verify that eigencube_gui renders frames with active solution and progress info."""
        import eigencube_gui

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        mock_solution = [((1, 0, 0), 1), ((0, 1, 0), -1)]
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "preview_solution.png")
            result = eigencube_gui.render_frame_to_image(
                solved_cube,
                out_path,
                solution=mock_solution,
                move_index=1,
                current_move=mock_solution[0],
            )
            self.assertEqual(result, out_path)
            self.assertTrue(os.path.exists(out_path))
            self.assertGreater(os.path.getsize(out_path), 0)

    def test_gui_window_icon_loads_without_sdl_image(self):
        """Icon must load even when pygame can only decode BMP (no SDL_image), via the Pillow fallback."""
        import pygame
        import eigencube_gui

        with Image.open(eigencube_gui.REPO_DIR / "img" / "icon.png") as source:
            w, h = source.size
            samples = ((0, 0), (w // 2, h // 2), (w - 1, h - 1))
            expected_alpha = "A" in source.getbands() or "transparency" in source.info  # mirrors the code's probe
            rgba_source = source.convert("RGBA")  # get_at always yields RGBA, whatever the source stores
            expected_pixels = [rgba_source.getpixel(p) for p in samples]
        with patch("pygame.image.load", side_effect=pygame.error("File is not a Windows BMP file")):
            icon = eigencube_gui.load_window_icon()
        self.assertEqual(icon.get_size(), (w, h))
        self.assertEqual(bool(icon.get_flags() & pygame.SRCALPHA), expected_alpha)
        for point, expected in zip(samples, expected_pixels):  # pixel bytes preserved
            self.assertEqual(tuple(icon.get_at(point)), expected)

    def test_gui_window_icon_fallback_keeps_opaque_rgb(self):
        """An alpha-free icon stays an alpha-free surface via the Pillow fallback, like a full pygame build."""
        import pygame
        import eigencube_gui

        with tempfile.TemporaryDirectory() as tmp_dir:
            img_dir = os.path.join(tmp_dir, "img")
            os.makedirs(img_dir)
            Image.new("RGB", (16, 8), (10, 20, 30)).save(os.path.join(img_dir, "icon.png"))
            with patch.object(eigencube_gui, "REPO_DIR", Path(tmp_dir)), \
                 patch("pygame.image.load", side_effect=pygame.error("File is not a Windows BMP file")):
                icon = eigencube_gui.load_window_icon()
        self.assertEqual(icon.get_size(), (16, 8))
        self.assertFalse(icon.get_flags() & pygame.SRCALPHA)
        self.assertEqual(tuple(icon.get_at((8, 4))), (10, 20, 30, 255))  # opaque pixels, read back as RGBA

    def test_gui_window_icon_fallback_keeps_palette_alpha(self):
        """A palette (P-mode) icon with a tRNS chunk must keep its transparency via the Pillow fallback."""
        import pygame
        import eigencube_gui

        with tempfile.TemporaryDirectory() as tmp_dir:
            img_dir = os.path.join(tmp_dir, "img")
            os.makedirs(img_dir)
            palette = Image.new("P", (8, 4))
            palette.putpalette([255, 0, 0, 0, 255, 0] + [0, 0, 0] * 254)
            palette.putpixel((1, 0), 1)  # the only opaque-colored-and-transparent index
            palette.info["transparency"] = 1
            palette.save(os.path.join(img_dir, "icon.png"))
            with patch.object(eigencube_gui, "REPO_DIR", Path(tmp_dir)), \
                 patch("pygame.image.load", side_effect=pygame.error("File is not a Windows BMP file")):
                icon = eigencube_gui.load_window_icon()
        self.assertTrue(icon.get_flags() & pygame.SRCALPHA)
        self.assertEqual(tuple(icon.get_at((0, 0))), (255, 0, 0, 255))
        self.assertEqual(tuple(icon.get_at((1, 0))), (0, 255, 0, 0))

    def test_gui_frame_saves_without_sdl_image(self):
        """Frames must save even when pygame can only decode BMP (no SDL_image), via the Pillow fallback."""
        import pygame
        import eigencube_gui

        surf = pygame.Surface((40, 20))
        surf.set_at((5, 5), (255, 0, 0))
        surf.set_at((30, 15), (0, 255, 0))
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "frame.png")
            with patch("pygame.image.save", side_effect=NotImplementedError("saving images of extended format is not available")):
                eigencube_gui.save_frame(surf, out_path)
            with Image.open(out_path) as saved:
                self.assertEqual(saved.mode, "RGB")  # opaque surface stays opaque, like pygame's own save
                self.assertEqual(saved.size, surf.get_size())  # dimensions preserved
                self.assertEqual(saved.getpixel((5, 5)), (255, 0, 0))  # pixels preserved
                self.assertEqual(saved.getpixel((30, 15)), (0, 255, 0))

    def test_gui_frame_saves_keep_alpha_surface_rgba(self):
        """An alpha surface must export as RGBA via the Pillow fallback, like pygame's own save."""
        import pygame
        import eigencube_gui

        surf = pygame.Surface((40, 20), pygame.SRCALPHA)
        surf.set_at((5, 5), (255, 0, 0, 128))
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "frame.png")
            with patch("pygame.image.save", side_effect=NotImplementedError("saving images of extended format is not available")):
                eigencube_gui.save_frame(surf, out_path)
            with Image.open(out_path) as saved:
                self.assertEqual(saved.mode, "RGBA")
                self.assertEqual(saved.getpixel((5, 5)), (255, 0, 0, 128))
                self.assertEqual(saved.getpixel((30, 15)), (0, 0, 0, 0))

    def test_gui_text_bubble_and_buttons(self):
        """Verify GUI button and text bubble components render without error across edge cases."""
        import eigencube_gui

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        eigencube_gui.init_display()

        # Button creation
        surf, rect = eigencube_gui.create_button("Test Button", 10, 20, 120, 35, (0, 0, 0), (255, 255, 255))
        self.assertEqual(surf.get_size(), (120, 35))
        self.assertEqual(rect.topleft, (10, 20))

        # Text bubble edge cases: empty text, bold prefix, progress bar
        h_plain = eigencube_gui.draw_text_bubble("Plain text", 10, 10, 200)
        h_bold = eigencube_gui.draw_text_bubble("Prefix: remaining text", 10, 10, 200, progress=0.5, bold_part="Prefix:")
        h_empty = eigencube_gui.draw_text_bubble("", 10, 10, 200, progress=1.0)
        self.assertGreater(h_plain, 0)
        self.assertGreater(h_bold, 0)
        self.assertGreater(h_empty, 0)

    def test_gui_font_fallback(self):
        """Verify font loading falls back safely to system fonts if font files cannot be loaded."""
        from unittest.mock import patch
        import eigencube_gui

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        orig_regular, orig_bold = eigencube_gui.font_regular, eigencube_gui.font_bold
        orig_max_h = eigencube_gui.MAX_TEXT_HEIGHT
        try:
            eigencube_gui.font_regular = None
            eigencube_gui.font_bold = None
            with patch("pygame.font.Font", side_effect=Exception("Simulated missing font")):
                eigencube_gui.init_display()
                self.assertIsNotNone(eigencube_gui.font_regular)
                self.assertIsNotNone(eigencube_gui.font_bold)
                surf, _ = eigencube_gui.create_button("Fallback", 0, 0, 100, 30, (0, 0, 0), (255, 255, 255))
                self.assertEqual(surf.get_size(), (100, 30))
        finally:
            eigencube_gui.font_regular, eigencube_gui.font_bold = orig_regular, orig_bold
            eigencube_gui.MAX_TEXT_HEIGHT = orig_max_h

    def test_gui_main_lifecycle(self):
        """Verify GUI main loop initializes, creates all UI components, and exits cleanly on QUIT."""
        import pygame
        import eigencube_gui

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        pygame.init()
        pygame.event.post(pygame.event.Event(pygame.QUIT))
        eigencube_gui.main()

    def test_gui_lifecycle_reinit_after_quit(self):
        """Verify re-initialization after pygame.quit() does not retain stale font handles or crash."""
        import pygame
        import eigencube_gui

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        eigencube_gui.init_display()
        pygame.quit()
        # Second session: re-init and ensure font rendering and button creation succeed without segfault
        eigencube_gui.init_display()
        surf, rect = eigencube_gui.create_button("Reinit Test", 0, 0, 100, 30, (0, 0, 0), (255, 255, 255))
        self.assertEqual(surf.get_size(), (100, 30))
        h = eigencube_gui.draw_text_bubble("Reinit Bubble", 0, 0, 200)
        self.assertGreater(h, 0)

    def test_gui_assets_independent_of_working_directory(self):
        """Verify the GUI starts from any directory and loads its bundled assets, not system fallbacks."""
        import pygame
        import eigencube_gui
        from unittest import mock

        os.environ["SDL_VIDEODRIVER"] = "dummy"
        pygame.quit()
        original_cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as elsewhere:
            try:
                os.chdir(elsewhere)
                with mock.patch("pygame.font.SysFont", side_effect=AssertionError("bundled font not found")):
                    eigencube_gui.init_display()
                out_path = eigencube_gui.render_frame_to_image(solved_cube, os.path.join(elsewhere, "frame.png"))
                self.assertTrue(os.path.exists(out_path))
            finally:
                os.chdir(original_cwd)

    def test_gui_offscreen_surface_preservation(self):
        """Verify passing a custom surface to init_display is preserved across subsequent draw calls."""
        import pygame
        import eigencube_gui

        custom_surf = pygame.Surface((eigencube_gui.WIDTH, eigencube_gui.HEIGHT))
        eigencube_gui.init_display(surface=custom_surf)
        self.assertIs(eigencube_gui.screen, custom_surf)
        eigencube_gui.draw_cube_static(solved_cube)
        self.assertIs(eigencube_gui.screen, custom_surf)
        eigencube_gui.draw_text_bubble("Test Offscreen", 10, 10, 200)
        self.assertIs(eigencube_gui.screen, custom_surf)


class TestLearnedSteps(unittest.TestCase):
    """Steps are learned at runtime, so their invariants must hold for any sequence."""

    # A commutator of two adjacent faces: disturbs few cubelets.
    commutator = (((-1, 0, 0), -1), ((0, -1, 0), -1), ((-1, 0, 0), 1), ((0, -1, 0), 1))

    def setUp(self):
        import eigencube
        self.steps = eigencube.steps
        self.saved_steps = dict(self.steps)

    def tearDown(self):
        self.steps.clear()
        self.steps.update(self.saved_steps)

    def test_symmetric_sequences_act_as_conjugates(self):
        """Under every cube symmetry Q, the image of a sequence acts as Q * sequence * Q^T."""
        from eigencube import symmetries, symmetric_move
        rng = np.random.default_rng(0)
        sequence = tuple(moves[i] for i in rng.integers(len(moves), size=8))
        original = dict(apply_step_to_cube(sequence, solved_cube))
        for Q in symmetries:
            image = dict(apply_step_to_cube(tuple(symmetric_move(Q, m) for m in sequence), solved_cube))
            for cubelet, rotation in original.items():
                image_cubelet = tuple(int(x) for x in Q @ cubelet)
                np.testing.assert_array_equal(image[image_cubelet], Q @ np.array(rotation) @ Q.T)

    def assert_step_table_invariants(self):
        from eigencube import num_disturbed, inverse_step, MAX_MACRO_DISTURBANCE
        single_moves = {apply_step_to_cube((m,), solved_cube) for m in moves}
        for effect, step in self.steps.items():
            # Keyed by its true effect, and within the disturbance cap.
            self.assertEqual(apply_step_to_cube(step, solved_cube), effect)
            self.assertTrue(0 < num_disturbed(effect) <= MAX_MACRO_DISTURBANCE)
            # Never a longer duplicate of a single move.
            self.assertTrue(len(step) == 1 or effect not in single_moves)
            # The inverse is a known step too.
            self.assertIn(apply_step_to_cube(inverse_step(step), solved_cube), self.steps)

    def test_step_table_starts_with_exactly_the_single_moves(self):
        self.assertEqual(sorted(self.saved_steps.values()), sorted((m,) for m in moves))

    def test_learned_steps_satisfy_table_invariants(self):
        """Learning keeps every entry keyed by its effect, capped, deduplicated and closed under inverse."""
        from eigencube import learn_step
        learn_step(self.commutator)
        self.assertIn(apply_step_to_cube(self.commutator, solved_cube), self.steps)
        # Search results that equal a single move must not duplicate it.
        for move in moves:
            learn_step((move,))
            learn_step(5 * (move,))
            learn_step((move, move))
        self.assert_step_table_invariants()

    def test_step_matches_move_by_move_application(self):
        """Applying a step equals applying its moves one by one, from any state."""
        from eigencube import learn_step, shuffle
        learn_step(self.commutator)
        cube = shuffle(solved_cube, iterations=50, seed=7)
        for step in self.steps.values():
            expected = cube
            for move in step:
                expected = apply_move_to_cube(move, expected)
            self.assertEqual(apply_step_to_cube(step, cube), expected)

    def test_learn_step_rejects_trivial_and_disruptive_sequences(self):
        """Sequences that disturb nothing, or too much, are not worth learning."""
        from eigencube import learn_step
        learn_step(())
        learn_step((moves[0], inverse_move(moves[0])))  # Disturbs nothing.
        learn_step(tuple(m for m in moves if m[1] == 1))  # Disturbs 20 of 26 cubelets.
        self.assertEqual(self.steps, self.saved_steps)

    def test_corner_twist_phase_solves_pure_twists(self):
        """The endgame solves a cube whose only defect is two twisted corners."""
        from eigencube import num_disturbed
        left, top, bottom = ((0, -1, 0), 1), ((0, 0, 1), 1), ((0, 0, -1), 1)
        twist_one_corner = 2 * (inverse_move(left), inverse_move(top), left, top)
        twist_two_corners = twist_one_corner + (bottom,) + twist_one_corner * 2 + (inverse_move(bottom),)
        cube = apply_step_to_cube(twist_two_corners, solved_cube)
        self.assertEqual(num_disturbed(cube), 2)
        solution = solve(cube)
        self.assertTrue(is_cube_solved(apply_step_to_cube(tuple(solution), cube)))


if __name__ == "__main__":
    unittest.main()
