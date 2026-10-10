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
    norm1,
    rotation_matrix,
    describe_move,
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

    def test_moves_turn_the_way_they_are_named(self):
        """Direction 1 is clockwise as seen looking at the face `v` points to, for every face.

        Seen from outside the face, a clockwise quarter turn about v takes each vector p
        perpendicular to v to M p with v · (p × M p) = -|p|²; counterclockwise gives +|p|².
        Order-4, det and axis checks all pass for a mirror-imaged turn, so only this catches it.
        """
        for move in moves:
            v, direction = move
            M = rotation_matrix(move)
            for p in unit_vectors:
                if np.dot(p, v) == 0:
                    self.assertEqual(np.dot(v, np.cross(p, M @ np.array(p))), -direction,
                                     f"{describe_move(move)} turns the wrong way")

    def test_opposite_face_turns_the_other_way(self):
        """Turning a face clockwise is the same rotation as turning its opposite face counterclockwise."""
        for (x, y, z), direction in moves:
            self.assertTrue(np.array_equal(rotation_matrix(((x, y, z), direction)),
                                           rotation_matrix(((-x, -y, -z), -direction))))

    def test_standard_notation(self):
        """Clockwise turns of U D F B R L move cubelets where Singmaster notation says they go."""
        clockwise = {  # face: (a cubelet's home, where the clockwise turn takes it)
            (0, 0, 1): ((1, 0, 1), (0, -1, 1)),     # U: front-top edge to left-top
            (0, 0, -1): ((1, 0, -1), (0, 1, -1)),   # D: front-bottom edge to right-bottom
            (1, 0, 0): ((1, 0, 1), (1, 1, 0)),      # F: top-front edge to right-front
            (-1, 0, 0): ((-1, 0, 1), (-1, -1, 0)),  # B: top-back edge to left-back
            (0, 1, 0): ((1, 1, 0), (0, 1, 1)),      # R: front-right edge to top-right
            (0, -1, 0): ((1, -1, 0), (0, -1, -1)),  # L: front-left edge to bottom-left
        }
        for v, (home, after) in clockwise.items():
            self.assertEqual(tuple(rotation_matrix((v, 1)) @ np.array(home)), after,
                             describe_move((v, 1)))

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
        """When random_weight=0, exceeding max_moves returns None instead of searching indefinitely."""
        # Scramble with 3 moves
        cube = solved_cube
        scramble = [((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 1)]
        for m in scramble:
            cube = apply_move_to_cube(m, cube)

        # Budget of 2 simulated moves is insufficient to solve a 3-move scramble
        res = astar(cube, is_cube_solved, apply_move_to_cube, random_weight=0, max_moves=2)
        self.assertIsNone(res)

        # Sufficient budget succeeds
        res = astar(cube, is_cube_solved, apply_move_to_cube, random_weight=0, max_moves=5000)
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

        # Initial budget of 5 moves will trigger restarts but budget expands by 1.5x until solved
        res = astar(cube, is_cube_solved, apply_move_to_cube, random_weight=0.25, max_moves=5)
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
            get_moves=lambda s: [],
            random_weight=0.25,
            max_moves=1000,
        )
        self.assertIsNone(res)

        # 2. Finite 3-state cyclic component {0, 1, 2} with unreachable goal 99
        # Moves simulated will reach max_moves=2, but once frontier empties it must return None
        dummy_move = ((1, 0, 0), 1)
        res_cyclic = astar(
            0,
            lambda x: x == 99,
            lambda m, s: (s + 1) % 3,
            get_moves=lambda s: [dummy_move],
            random_weight=0.25,
            max_moves=2,
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


if __name__ == "__main__":
    unittest.main()

