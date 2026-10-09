"""The Eigencube explainer, one Manim scene per chapter. Narration lives right next to its animation."""

from kit import *

FRONT, RIGHT_FACE, TOP = (1, 0, 0), (0, 1, 0), (0, 0, 1)
TURN_TOP, TURN_FRONT = (TOP, 1), (FRONT, 1)
CORNER, EDGE, CENTER = (1, 1, 1), (1, 0, 1), (0, 0, 1)
CAMERA = dict(phi=64 * DEGREES, theta=28 * DEGREES, zoom=1.2)
DEMO_MOVES = [TURN_TOP, TURN_FRONT, ((0, 1, 0), -1), ((0, 0, -1), 1)]


def sticker_arrows(cubelet, length=1.7 * SPACING):
    """One arrow per column of diag(c), from the cubelet's center out through its sticker."""
    return VGroup(*(piercing_arrow(at(cubelet), at(cubelet) + 0.5 * SPACING * n,
                                   at(cubelet) + length * n, STICKER[ints(n)])
                    for n in np.diag(cubelet).T if n.any()))


def assert_twisted_home(state):
    """What the film says of our corner after a top and a front turn: R c = c, yet it isn't
    solved, since R diag(c) != diag(c)."""
    rotation = dict(state)[CORNER]
    assert eigencube.position(CORNER, rotation) == CORNER
    assert not eigencube.is_cubelet_solved(CORNER, rotation)


def sticker_colors(cubelet):
    """The colors of sticker_arrows(cubelet), in order (to highlight each in its own color)."""
    return [STICKER[ints(n)] for n in np.diag(cubelet).T if n.any()]


def arm(direction, length, color, thickness=0.04):
    """An arrow from the core along `direction` (a unit address), out of the cube.

    Where it leaves through a face the camera sees (+x, +y, +z), the stretch outside is drawn
    over the cube with an outline; on the far side the cube rightly hides it. Each arm starts at
    the core's face, not its center: Cairo depth-sorts the stretch inside the core badly, and lets
    bits of it show through the cube."""
    start, tip = at(direction, 0.5), at(direction, length)
    if min(direction) < 0:
        return SolidArrow(start, tip, color, thickness=thickness)
    return piercing_arrow(start, at(direction, 1.55), tip, color, thickness=thickness)


def screen_right(theta):
    """The direction in the world that points right on screen, for a camera turned by `theta`."""
    return np.array([-np.sin(theta), np.cos(theta), 0])


def vector_2d(direction, color):
    """A flat arrow with a thin black outline, for the plane of the linear algebra review."""
    return Vector(direction, color=color).set_stroke(BLACK, width=11, background=True)


def position_arrow(point, scale=1.0, start=ORIGIN):
    """The position vector from the core to cubelet address `point`.

    It is drawn over the cube rather than depth-sorted into it: it runs through the cubelets, and
    sorted against them it would show only in fragments."""
    return overlay_arrow(start, at(point, scale), POSITION, thickness=0.04)


def vector_tex(v, **kwargs):
    return MathTex(vector_tex_string(v), **kwargs)


def vector_tex_string(v):
    return r"(%s)" % r",\, ".join(str(int(x)) for x in v)


def stacked_tex(rows, colors, **kwargs):
    """Rows of math typeset as one formula, aligned at their `&`, one submobject per row.

    Typeset together, the rows sit on evenly spaced baselines; separate formulas stacked by their
    bounding boxes would not be, since subscripts and descenders push some of them apart."""
    tex = MathTex(*(row + r" \\[4pt]" for row in rows[:-1]), rows[-1], **kwargs)
    for row, color in zip(tex, colors):
        row.set_color(color)
    return tex


def address_tex(c, **kwargs):
    """The label c = (x, y, z) for cubelet address `c`."""
    return MathTex("c = " + vector_tex_string(c), **kwargs)


def matrix_tex(rows, column_colors=None):
    """A bracketed integer matrix (or column vector) with optionally colored columns."""
    m = IntegerMatrix([[int(v) for v in row] for row in rows], h_buff=0.85)
    for column, color in zip(m.get_columns(), column_colors or []):
        column.set_color(color)
    return m


def color_name(vector):
    return eigencube.color_names[vector].lower()


def column_colors(cubelet):
    """Each column of diag(c) in its sticker's color (gray where there is no sticker)."""
    return [STICKER[ints(n)] if n.any() else GREY_D for n in np.diag(cubelet).T]


def rotation_rule(font_size):
    """The review's rule, which the film recalls word for word whenever it leans on it."""
    return VGroup(*(MathTex(line, font_size=font_size) for line in [
        r"\text{a rotation about the origin} \;=\; \text{a matrix } R",
        r"R\,\vec{u} \;=\; \vec{u}\,' \quad (\vec{u} \text{ rotated})",
    ])).arrange(DOWN, aligned_edge=LEFT, buff=0.3 * font_size / 40)


def labeled_matrix(label, rows, colors):
    return VGroup(MathTex(label, font_size=40), matrix_tex(rows, colors)).arrange(RIGHT)


# --- 1. The sticker nightmare -------------------------------------------------

def sticker_slots(cube):
    """Where each sticker currently is: (cubelet, color) -> (position R c, facing R n)."""
    return {(c, ints(n)): (ints(eigencube.position(c, r)), ints(np.array(r) @ n))
            for c, r in cube for n in np.diag(c).T if n.any()}


def sticker_list_layout():
    """Every sticker location of the cube as (position, direction), in a flat-list order."""
    face_order = [TOP, (-1, 0, 0), (0, -1, 0), FRONT, RIGHT_FACE, (0, 0, -1)]
    return sorted(sticker_slots(eigencube.solved_cube).values(),
                  key=lambda s: (face_order.index(s[1]), s[0]))


class StickerNightmare(Narrated):
    def construct(self):
        self.set_camera_orientation(**{**CAMERA, "zoom": 1.35})
        cube = CubeMobject(eigencube.shuffle(eigencube.solved_cube, 25, seed=5))
        self.begin_ambient_camera_rotation(rate=0.12)
        with self.voice("Here's a Rubik's cube. Twenty-six little plastic cubelets, fifty-four "
                        "colored stickers, and about forty-three quintillion ways to mix them up."):
            # On screen from the first frame, which is what a player shows before it starts.
            self.add(cube)
            for move in [(TOP, 1), ((0, 1, 0), -1), (FRONT, 1), ((0, 0, -1), 1)]:
                self.play(cube.turn(move), run_time=0.7)
        self.stop_ambient_camera_rotation()

        with self.voice("Suppose you want to teach a computer this puzzle. The obvious first idea "
                        "is to write down the stickers: fifty-four colors in a list."):
            layout = sticker_list_layout()
            width = 0.235
            cells = VGroup(*(Square(width * 0.9, fill_color=STICKER[d], fill_opacity=1,
                                    stroke_width=0) for _, d in layout))
            cells.arrange(RIGHT, buff=width * 0.1).move_to(DOWN * 0.3)
            indices = VGroup(*(Text(str(i), font_size=12, color=GREY_B).next_to(cell, DOWN, 0.08)
                               for i, cell in enumerate(cells) if i % 9 == 0))
            brackets = VGroup(Text("[", font_size=40).next_to(cells, LEFT, 0.05),
                              Text("]", font_size=40).next_to(cells, RIGHT, 0.05))
            listing = self.hud(VGroup(brackets, cells, indices))
            self.play(FadeOut(cube), run_time=0.8)
            self.play(LaggedStart(*(FadeIn(c, shift=0.2 * UP) for c in cells), lag_ratio=0.03),
                      FadeIn(brackets), FadeIn(indices), run_time=2.5)

        with self.voice("And then the pain begins. A single quarter turn of the top face sends "
                        "twenty of those stickers to new slots, scattered all over the list."):
            slot = {s: i for i, s in enumerate(layout)}
            before = sticker_slots(eigencube.solved_cube)
            after = sticker_slots(eigencube.apply_move_to_cube(TURN_TOP, eigencube.solved_cube))
            moved = [(cells[slot[before[s]]], cells[slot[after[s]]].get_center())
                     for s in before if before[s] != after[s]]
            assert len(moved) == 20
            arcs = VGroup(*(ArcBetweenPoints(cell.get_center() + 0.15 * UP, target + 0.15 * UP,
                                             angle=-PI / 2.2, color=ACCENT, stroke_width=1.5)
                            for cell, target in moved))
            self.hud(arcs, overlay=True)
            self.wait(1.2)
            self.play(LaggedStart(*(Draw(a) for a in arcs), lag_ratio=0.06), run_time=2.5)
            self.play(*(cell.animate(path_arc=-PI / 2).move_to(target) for cell, target in moved),
                      run_time=2)
            self.play(FadeOut(arcs))

        with self.voice("Twelve possible moves means twelve permutation tables. Then you need more "
                        "bookkeeping to know which stickers share a corner, how corners twist, "
                        "and how edges flip. It works, but the geometry is gone. It's just index "
                        "shuffling."):
            tables = VGroup(*(Text(f"move_{i:02d} = [ ... ]", font="DejaVu Sans Mono",
                                   font_size=20, color=GREY_B) for i in range(12)))
            tables.arrange_in_grid(4, 3, buff=(0.6, 0.25)).to_edge(UP, buff=0.6)
            self.hud(tables)
            self.play(LaggedStart(*(FadeIn(t) for t in tables), lag_ratio=0.15), run_time=3)
            extra = self.hud(Text("+ which stickers share a corner, corner twists, edge flips",
                                  font_size=24, color=GREY_A).next_to(listing, DOWN, 0.7))
            self.when_said("bookkeeping")
            self.play(Write(extra))

        logo = SVGMobject(str(HERE.parent / "img" / "logo.svg")).scale(1.1)
        title = Text("Eigencube", font_size=80, weight=BOLD)
        subtitle = Text("a Rubik's cube made of linear algebra", font_size=34, color=ACCENT)
        title_card = self.hud(VGroup(logo, title, subtitle).arrange(DOWN, buff=0.3).shift(0.8 * UP))
        with self.voice("Today I want to show you a different way. A way where the entire cube, its "
                        "state, its moves, even the test for whether it's solved, is nothing but "
                        "vectors, matrices, and dot products."):
            self.play(FadeOut(listing), FadeOut(tables), FadeOut(extra), run_time=1)
            self.play(FadeIn(logo, scale=0.8), Write(title), run_time=1.5)
            self.play(FadeIn(subtitle, shift=0.2 * UP))

        with self.voice("It's called Eigencube. The whole model, plus a solver, fits in under four "
                        "hundred lines of Python, with no lookup tables at all. Let's see how it "
                        "works."):
            facts = self.hud(VGroup(
                Text("eigencube.py  ·  < 400 lines  ·  no lookup tables",
                     font="DejaVu Sans Mono", font_size=24, color=GREY_B),
                Text("github.com/smolkaj/eigencube", font="DejaVu Sans Mono", font_size=24,
                     color=ACCENT),
            ).arrange(DOWN, buff=0.2).next_to(title_card, DOWN, 0.4))
            self.play(FadeIn(facts, shift=0.2 * UP))
        self.play(FadeOut(title_card), FadeOut(facts))


# --- 2. A short review of linear algebra ------------------------------------------

class LinearAlgebraReview(Narrated):
    def construct(self):
        name = Text("Essence of Linear Algebra", font_size=54, weight=BOLD)
        by = Text("3Blue1Brown  ·  Grant Sanderson", font_size=30, color=ACCENT)
        card = self.label(VGroup(name, by).arrange(DOWN, buff=0.3))
        with self.voice("Before we dive in, a shoutout. This whole project grew out of a video "
                        "series by 3Blue1Brown: Essence of Linear Algebra, by Grant Sanderson. "
                        "We'll lean on a few of its ideas, so let's take a minute to review them."):
            self.play(Write(name), run_time=1.5)
            self.play(FadeIn(by, shift=0.2 * UP))
        card.add_background_rectangle(opacity=0.9, buff=0.2)  # Keeps it legible over the grid.
        self.play(card.animate.scale(0.5).to_corner(UL))

        plane = NumberPlane(x_range=[-10, 10], y_range=[-10, 10],
                            background_line_style={"stroke_color": BLUE_D, "stroke_opacity": 0.6})
        v = vector_2d([2, 1], POSITION)
        v_label = MathTex(r"\begin{bmatrix} 2 \\ 1 \end{bmatrix}", color=POSITION).next_to(
            v.get_end(), RIGHT, 0.15)
        self.label(v_label)
        with self.voice("Here's a flat plane with a grid. A vector is an arrow that starts at the "
                        "origin, and we describe it by its coordinates. This one is two, one: two "
                        "steps to the right, and one step up."):
            self.play(Draw(plane, lag_ratio=0.05), run_time=2)
            self.play(GrowArrow(v))
            self.play(Write(v_label))

        e_x, e_y = vector_2d(RIGHT, X_COLOR), vector_2d(UP, Y_COLOR)
        x_label = MathTex(r"\mathbf{e}_x", color=X_COLOR).next_to(e_x, DOWN, 0.15)
        y_label = MathTex(r"\mathbf{e}_y", color=Y_COLOR).next_to(e_y, LEFT, 0.15)
        self.label(x_label, y_label)
        with self.voice("Those steps are measured with two special vectors. e x is one step to "
                        "the right, drawn here in green. e y is one step up, drawn in red. "
                        "Together, they're called the basis. Essence of Linear Algebra calls them "
                        "i-hat and j-hat."):
            self.when_said("e x is")
            self.play(GrowArrow(e_x), Write(x_label))
            self.when_said("e y is")
            self.play(GrowArrow(e_y), Write(y_label))

        steps = VGroup(vector_2d(RIGHT, X_COLOR).shift(RIGHT),
                       vector_2d(UP, Y_COLOR).shift(2 * RIGHT))
        combination = MathTex(r"\vec{u}", "=", r"2\,\mathbf{e}_x", "+", r"1\,\mathbf{e}_y",
                              font_size=48).to_corner(UR).add_background_rectangle()
        combination[3].set_color(X_COLOR)  # Indices count the background rectangle first.
        combination[5].set_color(Y_COLOR)
        self.label(combination)
        with self.voice("Every vector is a combination of them. Ours is two steps of e x, plus one "
                        "step of e y."):
            self.when_said("two steps")
            self.play(TransformFromCopy(e_x, steps[0]))
            self.when_said("one step")
            self.play(TransformFromCopy(e_y, steps[1]))
            self.play(Write(combination))
        self.play(FadeOut(steps))

        quarter_turn = [[0, -1], [1, 0]]
        with self.voice("Now, a matrix is a way to transform all of space. Let's rotate everything "
                        "a quarter turn, counterclockwise, about the origin."):
            self.play(FadeOut(x_label), FadeOut(y_label), FadeOut(v_label))
            self.when_said("rotate")
            self.play(Rotate(VGroup(plane, e_x, e_y, v), angle=PI / 2, about_point=ORIGIN),
                      run_time=3)

        x_lands = MathTex(r"\begin{bmatrix} 0 \\ 1 \end{bmatrix}", color=X_COLOR).next_to(
            e_x.get_end(), RIGHT, 0.15).add_background_rectangle()
        y_lands = MathTex(r"\begin{bmatrix} -1 \\ 0 \end{bmatrix}", color=Y_COLOR).next_to(
            e_y.get_end(), DOWN, 0.15).add_background_rectangle()
        m = IntegerMatrix(quarter_turn, h_buff=0.9).add_background_rectangle()
        m.get_columns()[0].set_color(X_COLOR)
        m.get_columns()[1].set_color(Y_COLOR)
        m_group = VGroup(MathTex("R =", font_size=48), m).arrange(RIGHT).next_to(
            combination, DOWN, 0.5, aligned_edge=RIGHT)
        self.label(x_lands, y_lands, m_group)
        with self.voice("Here's the beautiful part. e x has landed on zero, one. e y has landed on "
                        "minus one, zero. Write those landing spots side by side, as columns, and "
                        "that is the matrix of this rotation."):
            self.when_said("landed on zero")
            self.play(Write(x_lands))
            self.when_said("landed on minus")
            self.play(Write(y_lands))
            self.when_said("side by side")
            self.play(Write(m_group[0]), FadeIn(m.get_brackets()), FadeIn(m.background_rectangle))
            # Plain fades: a transform from a copy would leave the copy, not the matrix, on screen.
            self.play(Indicate(x_lands, color=X_COLOR), FadeIn(m.get_columns()[0], shift=0.4 * RIGHT))
            self.play(Indicate(y_lands, color=Y_COLOR), FadeIn(m.get_columns()[1], shift=0.4 * RIGHT))

        product = MathTex(r"R \begin{bmatrix} 2 \\ 1 \end{bmatrix}", "=",
                          r"2 \begin{bmatrix} 0 \\ 1 \end{bmatrix}", "+",
                          r"1 \begin{bmatrix} -1 \\ 0 \end{bmatrix}", "=",
                          r"\begin{bmatrix} -1 \\ 2 \end{bmatrix}", font_size=44)
        product[2].set_color(X_COLOR)
        product[4].set_color(Y_COLOR)
        product[6].set_color(POSITION)
        product.add_background_rectangle().to_edge(DOWN, buff=1.3)
        self.label(product)
        with self.voice("And our vector? It's still two steps of e x plus one of e y, just using the "
                        "moved ones. That's exactly what multiplying the matrix by the vector "
                        "computes: the rotated vector, minus one, two."):
            self.when_said("two steps")
            self.play(Write(product[:6]), run_time=2.5)
            self.when_said("the rotated vector")
            self.play(Write(product[6:]), Indicate(v, color=POSITION))

        rule = rotation_rule(font_size=34)
        box = SurroundingRectangle(rule, color=ACCENT, buff=0.3).set_fill(BLACK, opacity=0.9)
        # Below the x axis, so that the rotated e_y, which ends up pointing left, stays in view.
        rule_card = self.label(VGroup(box, rule).move_to(1.0 * DOWN))
        with self.voice("So here's our rule for the rest of the video, now in three dimensions. A "
                        "rotation about the origin is a matrix, and multiplying a vector by that "
                        "matrix gives you the rotated vector. If you want to see why that works, "
                        "Essence of Linear Algebra is the place to go."):
            self.play(FadeOut(product), FadeOut(x_lands), FadeOut(y_lands))
            self.when_said("A rotation")
            self.play(FadeIn(box), Write(rule), run_time=2)  # Shown as rule_card.

        eigen = MathTex(r"\text{eigenvector, eigenvalue 1:}\;\; \text{a vector the transformation "
                        r"leaves in place}",
                        font_size=36, color=ACCENT)
        # Boxed like the rule card, and kept above the grid as it turns: Cairo otherwise draws
        # whatever is being animated on top.
        eigen_box = SurroundingRectangle(eigen, color=ACCENT, buff=0.15).set_fill(BLACK, 0.9)
        self.label(VGroup(eigen_box, eigen).next_to(rule_card, DOWN, 0.2).set_z_index(1))
        with self.voice("One last word we'll need. If a transformation leaves a vector exactly where "
                        "it was, that vector is called an eigenvector, with eigenvalue one. A "
                        "quarter turn of the plane leaves no arrow in place. But in three "
                        "dimensions, every rotation has one: its axis. Hold on to that thought."):
            self.when_said("eigenvector")
            self.play(FadeIn(eigen_box), Write(eigen), run_time=2)
            self.when_said("A quarter turn")  # Every arrow moves.
            self.play(FadeOut(rule_card), run_time=0.6)
            turning = VGroup(plane, e_x, e_y, v)
            self.play(Rotate(turning, angle=PI / 2, about_point=ORIGIN), run_time=2.5)
            # In three dimensions: the plane tilts into view, and turns about an axis that stays.
            self.when_said("But in three")
            keep = {turning, *turning.get_family()}
            self.play(*(FadeOut(m) for m in self.mobjects if m not in keep), run_time=0.5)
            # The flat arrows' black outlines would cut gaps into the grid's axes in perspective.
            for flat in (e_x, e_y, v):
                flat.set_stroke(width=0, background=True)
            # Close in, so that the arrows read and the plane's edge stays out of the picture.
            self.move_camera(phi=60 * DEGREES, theta=-60 * DEGREES, zoom=1.6, run_time=1.5)
            axis = SolidArrow(ORIGIN, 2 * OUT, ACCENT, thickness=0.04)
            self.when_said("its axis")
            self.play(Draw(axis), run_time=0.8)
            # A half turn: the vector ends pointing sideways, opposite to where it started, where a
            # quarter turn would leave it pointing down the screen, as if along the axis.
            self.play(Rotate(turning, angle=PI, about_point=ORIGIN), run_time=2.5)
        self.wait(0.5)
        self.play(*(FadeOut(m) for m in self.mobjects))


# --- 3. The fixed frame ---------------------------------------------------------

class FixedFrame(Narrated):
    def construct(self):
        self.set_camera_orientation(**CAMERA)
        cube = CubeMobject()
        core = cubelet_mobject((0, 0, 0))
        cross = [c for c in cube.pieces if eigencube.norm1(c) == 1]
        self.add(cube)
        with self.voice("Let's start with the physical cube. Take it apart, and you find its "
                        "skeleton: a solid three-dimensional cross. A core, with six center pieces "
                        "fixed to it."):
            self.begin_ambient_camera_rotation(rate=0.1)
            self.when_said("skeleton")
            self.add(core)
            self.play(cube.animate.ghosted(cross), run_time=2)
        self.stop_ambient_camera_rotation()
        self.move_camera(**CAMERA, run_time=1.5)

        with self.voice("Here's the key observation. Every move turns an outer face around one arm "
                        "of that cross. Watch the centers: they spin in place, but they never go "
                        "anywhere. The cross is fixed. Forever."):
            self.play(cube.animate.ghosted(cross, 0.3), run_time=0.6)
            made = self.turn_while_speaking(cube, DEMO_MOVES)
        self.undo(cube, made, run_time=0.15)  # Quickly: nothing is said meanwhile.

        # Both arms of the vertical axis, and the turns around them: top and bottom.
        axis = VGroup(arm(TOP, 3.6, ACCENT), arm((0, 0, -1), 3.4, ACCENT))
        # In words: the symbols for a face turn (M) and its axis (v) only arrive with the moves.
        eigen = self.hud(VGroup(*(Tex(line, font_size=40, color=ACCENT) for line in [
            "an eigenvector of", "every turn around it"])).arrange(DOWN, aligned_edge=RIGHT)
                         .to_corner(UR))
        with self.voice("In the language of our review: each arm of the cross stays exactly where it "
                        "is, under every turn around it, like these turns of the top and the "
                        "bottom. It's an eigenvector of those turns. And that is where the name "
                        "Eigencube comes from."):
            self.play(Draw(axis))
            made = [TURN_TOP, ((0, 0, -1), 1)]
            self.when_said("like these turns")
            for move in made:
                self.play(cube.turn(move))
            self.when_said("eigenvector")
            self.play(Write(eigen))
            made += self.turn_while_speaking(cube, made)
        self.undo(cube, made)
        self.play(FadeOut(axis), FadeOut(eigen))

        arrows, labels = self.axes()
        self.say("Something that never moves is exactly what we want for a coordinate system. Put "
                 "the origin at the core, and let the arms of the cross be the axes.")
        for axis_arrow, label, words in zip(arrows, labels, ["x points to the front,",
                                                             "y to the right,",
                                                             "and z to the top."]):
            with self.voice(words, pause=0.1):
                self.play(Draw(axis_arrow), FadeIn(label), run_time=1.1)
        self.say("Each axis runs both ways, out through the opposite center, too.")

        corner = (1, -1, 1)
        path = [overlay_arrow(a, b, color, thickness=0.04)  # Drawn over the cube.
                for a, b, color in [(ORIGIN, at(FRONT), X_COLOR),
                                    (at(FRONT), at((1, -1, 0)), Y_COLOR),
                                    (at((1, -1, 0)), at(corner), Z_COLOR)]]
        c_arrow = position_arrow(corner)
        c_label = self.hud(address_tex(corner, color=POSITION, font_size=44)
                           .to_corner(UL))
        with self.voice("Now every cubelet gets an address: the vector c, from the core to the "
                        "cubelet's center. Take this corner. To reach it, go one step to the "
                        "front, one step to the left, and one step up. So c equals one, minus "
                        "one, one."):
            self.play(cube.animate.ghosted(opacity=0.08), run_time=1)
            # See-through, so the arrow to its center shows.
            self.play(cube.pieces[corner].animate.set_opacity(0.45))
            for step, words in zip(path, ["to the front", "to the left", "one step up"]):
                self.when_said(words)
                self.play(Draw(step), run_time=0.8)
            self.when_said("So c equals")
            self.play(Draw(c_arrow), FadeIn(c_label))
        self.play(FadeOut(*path))

        # The right center and the top-right edge are both in plain view; labels sit clear of axes.
        others = [((0, 1, 0), (0, 1.6, -0.75)), ((0, 1, 1), (0, 0.9, 1.85))]
        other_arrows = [position_arrow(c) for c, _ in others]
        other_labels = [self.facing_camera(vector_tex(c, color=POSITION, font_size=36)
                                           .add_background_rectangle(opacity=0.75)
                                           .move_to(at(spot))) for c, spot in others]
        domain = self.hud(MathTex(r"c \in \{-1, 0, 1\}^3", font_size=44)
                          .next_to(c_label, DOWN, 0.4, aligned_edge=LEFT))
        with self.voice("The coordinates of every cubelet are minus one, zero, or one. A center, "
                        "like this one, sits at zero, one, zero. And an edge, like this one, at "
                        "zero, one, one."):
            self.play(Write(domain))
            for (c, _), a, label, words in zip(others, other_arrows, other_labels,
                                              ["A center", "And an edge"]):
                self.when_said(words)
                self.play(cube.pieces[c].animate.set_opacity(0.45), Draw(a), FadeIn(label))
        self.wait(2)


# --- 4. Counting stickers with the 1-norm ---------------------------------------

class CountingStickers(Narrated):
    def construct(self):
        self.set_camera_orientation(**CAMERA)
        axes = self.show_axes()
        cube = CubeMobject()
        core = cubelet_mobject((0, 0, 0))
        self.add(cube, core)
        slab_labels = {x: self.hud(MathTex(f"x = {x}", r"\text{: %s slice}" % where, color=X_COLOR,
                                           font_size=44).to_corner(UL))
                       for x, where in [(1, "front"), (0, "middle"), (-1, "back"), ("\\pm 1", "outer")]}
        with self.voice("Is this a good coordinate system? Here's an early sign. Look along the x "
                        "axis. The cube splits into three slices: x equals one, at the front; x "
                        "equals zero, in the middle; and x equals minus one, at the back."):
            shown = None
            for x, words in [(1, "x equals one"), (0, "x equals zero"), (-1, "x equals minus")]:
                self.when_said(words)
                self.play(cube.animate.ghosted([c for c in cube.pieces if c[0] == x], 0.1),
                          *([FadeOut(shown)] if shown else []), FadeIn(slab_labels[x]), run_time=1.2)
                shown = slab_labels[x]

        with self.voice("Cubelets in the middle slice, with x equal to zero, are inside along x: "
                        "nothing of theirs faces front or back."):
            self.play(cube.animate.ghosted([c for c in cube.pieces if c[0] == 0], 0.1),
                      FadeOut(shown), FadeIn(slab_labels[0]))
            shown = slab_labels[0]
        with self.voice("Cubelets in the outer slices sit on the surface. Each one shows a sticker "
                        "on that side."):
            self.play(cube.animate.ghosted([c for c in cube.pieces if c[0] != 0], 0.1),
                      FadeOut(shown), FadeIn(slab_labels["\\pm 1"]))
            shown = slab_labels["\\pm 1"]
        switch = self.hud(VGroup(
            MathTex(r"|x| = 1:", r"\;\text{a sticker facing front or back}", font_size=34),
            MathTex(r"|x| = 0:", r"\;\text{no sticker along } x", font_size=34),
        ).arrange(DOWN, aligned_edge=LEFT).to_corner(UL))
        with self.voice("So the absolute value of x is a little on-off switch: does this cubelet "
                        "show a sticker along the x axis? And the same goes for y, and for z."):
            self.play(FadeOut(shown))
            self.play(Write(switch), run_time=2)
            self.when_said("the same goes")
            self.play(cube.animate.ghosted(cube.pieces))

        norm = self.hud(VGroup(MathTex(r"\|c\|_1 = |x| + |y| + |z|", font_size=44),
                               MathTex(r"= \#\,\text{stickers}", font_size=44))
                        .arrange(DOWN, aligned_edge=LEFT).to_corner(UL))
        with self.voice("Add up the three switches, and you get what's called the Manhattan norm of "
                        "c. And it counts the stickers."):
            self.play(FadeOut(switch))
            self.when_said("Manhattan")
            self.play(Write(norm[0]))
            self.when_said("counts")
            self.play(Write(norm[1]))

        pieces = [core] + list(cube.pieces.values())
        homes = [(0, 0, 0)] + list(cube.pieces)
        # The axes step aside while the cube is taken apart; they would run through the pieces.
        self.play(*(p.animate.shift(at(c, 0.45)) for p, c in zip(pieces, homes)),
                  self.camera.zoom_tracker.animate.set_value(0.68), FadeOut(axes[0]),
                  FadeOut(*axes[1]), run_time=2)

        examples = [(0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)]
        names = ["core", "centers", "edges", "corners"]
        lines = ["Zero for the hidden core, at zero, zero, zero.",
                 "One for each of the six centers, like zero, zero, one.",
                 "Two for the twelve edges, like one, zero, one.",
                 "And three for the eight corners, like one, one, one."]
        rows = VGroup()
        for stickers, (c, name) in enumerate(zip(examples, names)):
            count = sum(eigencube.norm1(v) == stickers for v in eigencube.vectors)
            rows.add(MathTex(vector_tex_string(c), r":\;",
                             "+".join(str(abs(x)) for x in c), f"= {stickers}",
                             r"\;\to\;" + rf"\text{{{count} {name}}}", font_size=28))
        rows.arrange(DOWN, aligned_edge=LEFT, buff=0.25).next_to(norm, DOWN, 0.6, aligned_edge=LEFT)
        self.hud(*rows)
        example_arrow = None
        for stickers, (c, row, line) in enumerate(zip(examples, rows, lines)):
            highlight = [p for c2, p in cube.pieces.items() if eigencube.norm1(c2) == stickers]
            new_arrow = position_arrow(c, 1.45) if any(c) else None
            with self.voice(line, pause=0.3):
                # The core is a plain dark block with no stickers; tint it while it is the one named.
                core_look = core.animate.set_fill(ACCENT if stickers == 0 else BODY,
                                                  opacity=1 if stickers == 0 else 0.1)
                self.play(core_look, *(p.animate.set_opacity(1 if p in highlight else 0.1)
                                       for p in cube.pieces.values()),
                          *([FadeOut(example_arrow)] if example_arrow else []),
                          *([Draw(new_arrow)] if new_arrow else []),
                          FadeIn(row, shift=0.2 * RIGHT), run_time=1.2)
            example_arrow = new_arrow

        total = self.hud(MathTex(r"1 + 6 + 12 + 8 = 27", font_size=40, color=ACCENT)
                         .next_to(rows, DOWN, 0.4, aligned_edge=LEFT))
        with self.voice("One plus six plus twelve plus eight: all twenty-seven pieces, sorted by "
                        "type, with nothing more than a sum."):
            self.play(*(p.animate.set_opacity(1) for p in pieces), FadeOut(example_arrow),
                      Write(total))
        self.play(FadeOut(rows), FadeOut(total), FadeOut(norm))
        self.play(*(p.animate.shift(-at(c, 0.45)) for p, c in zip(pieces, homes)),
                  self.camera.zoom_tracker.animate.set_value(CAMERA["zoom"]), FadeIn(axes[0]),
                  FadeIn(*axes[1]), run_time=2)
        self.wait(0.5)


# --- 5. Colors are basis vectors ------------------------------------------------

class ColorsAreVectors(Narrated):
    def construct(self):
        self.set_camera_orientation(**CAMERA)
        cube = CubeMobject()  # No axes: the color arrows drawn below run along them.
        self.add(cube)
        centers = [c for c in cube.pieces if eigencube.norm1(c) == 1]
        with self.voice("Next question: how do we represent colors? A program might use strings, "
                        "or an enum. But look at the centers again while the cube turns. They "
                        "never move. And each one decides the color its whole face must have, "
                        "once the cube is solved."):
            made = self.turn_while_speaking(cube, DEMO_MOVES)

        # Long arrows, so that even the ones pointing away from the camera stand out; those longest,
        # since the cube hides most of their length. Orange's would reach into the address list.
        color_arrows = {c: arm(c, 3.6 if c in [(-1, 0, 0), (0, 0, -1)] else 3.0, STICKER[c],
                               thickness=0.045) for c in centers}
        # Each center's address, listed on screen as it is named: in the picture, the labels
        # would have to sit on the arrows and stickers they describe.
        order = [FRONT, RIGHT_FACE, TOP, (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
        addresses = self.hud(stacked_tex([r"&\text{%s center: }%s" % (color_name(c),
                                                                    vector_tex_string(c))
                                          for c in order], [STICKER[c] for c in order],
                                         font_size=32).to_corner(UL))
        tips = dict(zip(order, addresses))
        first = True
        for c, words in [(FRONT, "The green center always sits at one, zero, zero: one step "
                                 "along x."),
                         (RIGHT_FACE, "The red one, at zero, one, zero."),
                         (TOP, "And the white one, at zero, zero, one.")]:
            with self.voice(words):
                if first:  # Back to solved, under the first line.
                    self.undo(cube, made, run_time=0.2)
                    self.play(cube.animate.ghosted(centers, 0.12), run_time=0.8)
                    first = False
                self.play(Draw(color_arrows[c]), FadeIn(tips[c]), run_time=1.2)

        basis = self.hud(stacked_tex([r"\mathbf{e}_%s &= %s" % (axis, vector_tex_string(unit(i)))
                                      for i, axis in enumerate("xyz")], BASIS_COLORS,
                                     font_size=40).to_corner(UR))
        with self.voice("These three vectors have a name: the standard basis vectors, e x, e y, "
                        "and e z. One step along each axis: the same basis as in our review, plus "
                        "one more for the third dimension."):
            for b, words in zip(basis, ["e x", "e y", "and e z"]):
                self.when_said(words)
                self.play(Write(b), run_time=0.9)

        vectors = [unit(i) for i in range(3)] + [tuple(-x for x in unit(i)) for i in range(3)]
        palette = self.hud(stacked_tex(
            [r"%s\mathbf{e}_%s &= \text{%s}" % ("-" if any(x < 0 for x in v) else "", "xyz"[i % 3],
                                                color_name(v))
             for i, v in enumerate(vectors)], [STICKER[v] for v in vectors], font_size=40)
            .to_corner(UR))
        negatives = [v for v in vectors if any(x < 0 for x in v)]
        with self.voice("So why invent names for colors at all? A color simply is a unit vector. "
                        "Green is e x, red is e y, white is e z. And blue, orange, and yellow are "
                        "their opposites: minus e x, minus e y, and minus e z."):
            self.when_said("Green is")
            self.play(ReplacementTransform(basis, palette[:3]), run_time=1.5)
            self.when_said("And blue")
            self.play(*(Draw(color_arrows[c]) for c in negatives),
                      *(FadeIn(tips[c]) for c in negatives),
                      LaggedStart(*(Write(p) for p in palette[3:]), lag_ratio=0.4), run_time=2.5)
        # The opposites point away from the camera; look at them from the other side.
        # From behind, the arrows reach toward the address list, which steps aside meanwhile.
        with self.voice("They face away from us, so let's look from behind. There they are: blue, "
                        "orange, and yellow, each pointing the opposite way of its partner."):
            self.play(FadeOut(addresses))
            self.move_camera(phi=70 * DEGREES, theta=CAMERA["theta"] + PI, run_time=2.5)
        self.move_camera(**CAMERA, run_time=2)
        self.play(FadeIn(addresses))
        self.wait(1)


# --- 6. diag(c) -----------------------------------------------------------------

class DiagTrick(Narrated):
    def construct(self):
        # A wider angle than usual, so the arrow out of the front face does not point at the camera.
        self.set_camera_orientation(phi=66 * DEGREES, theta=34 * DEGREES, zoom=1.15)
        axes = self.show_axes()
        cube = CubeMobject()
        self.add(cube)
        c_label = self.hud(address_tex(CORNER, color=POSITION, font_size=44)
                           .to_corner(UL))
        with self.voice("Now for my favorite trick. Pick a cubelet: say, the front, right, top "
                        "corner. Its address is c equals one, one, one."):
            self.play(cube.animate.ghosted([CORNER], 0.1), run_time=1.5)
            self.when_said("Its address")
            self.play(FadeIn(c_label))

        decomposition = self.hud(MathTex(r"c", "=", r"1\,\mathbf{e}_x", "+", r"1\,\mathbf{e}_y",
                                         "+", r"1\,\mathbf{e}_z", font_size=44).to_corner(UL))
        for part, color in zip(decomposition[2::2], BASIS_COLORS):
            part.set_color(color)
        arrows = sticker_arrows(CORNER)
        with self.voice("Split that address along the axes: c is one e x, plus one e y, plus one "
                        "e z. Now draw those three pieces as arrows, starting from the cubelet "
                        "itself."):
            self.when_said("c is one")
            self.play(FadeOut(c_label))
            self.play(Write(decomposition), run_time=2)
            self.when_said("draw those")
            for a in arrows:
                self.play(Draw(a), run_time=0.8)

        self.say("Look where they go! Each one points straight out through one of the cubelet's "
                 "stickers. And since colors are basis vectors, each arrow is also the color of "
                 "the sticker it goes through.")
        for a, color, words in zip(arrows, sticker_colors(CORNER),
                                   ["Green, front.", "Red, right.", "White, top."]):
            with self.voice(words, pause=0.25):
                self.play(Indicate(a, color=color, scale_factor=1.4), run_time=1)

        diag = MathTex(r"\mathrm{diag}(c) =", font_size=44)
        matrix = matrix_tex(np.diag(CORNER), column_colors(CORNER))
        diag_group = self.hud(VGroup(diag, matrix).arrange(RIGHT).to_corner(UR))
        with self.voice("Now stack those three arrows side by side, as the columns of a matrix. "
                        "That's the diagonal matrix with c along its diagonal. We call it diag of "
                        "c."):
            self.when_said("side by side")
            self.play(Write(diag), FadeIn(matrix.get_brackets()))
            for column in matrix.get_columns():
                self.play(FadeIn(column, shift=0.3 * LEFT), run_time=0.7)

        # Each row is typeset as one formula, so that its words share a baseline whatever their
        # descenders (the subscript y would otherwise push "red" lower).
        symbols = MathTex(*(r"\mathbf{e}_%s \quad" % axis for axis in "xyz"), font_size=30)
        words = MathTex(*(r"\text{%s} \quad" % color_name(unit(i)) for i in range(3)),
                        font_size=30)
        VGroup(symbols, words).arrange(DOWN, buff=0.1).next_to(matrix, DOWN, 0.15)
        names = self.hud(VGroup(*(VGroup(symbol.set_x(column.get_x()), word.set_x(column.get_x()))
                                  .set_color(color)
                                  for symbol, word, column, color in zip(
                                      symbols, words, matrix.get_columns(), BASIS_COLORS))))
        self.say("Read its columns.")
        for column, name, words in zip(matrix.get_columns(), names,
                                       ["The first is one, zero, zero: that's e x, so green.",
                                        "The second is e y: red.",
                                        "And the third is e z: white."]):
            with self.voice(words, pause=0.2):
                self.play(Indicate(column, color=column.get_color()), FadeIn(name, shift=0.2 * DOWN))
        with self.voice("Each column is one sticker: the direction it faces, which is also its "
                        "color."):
            self.play(*(Indicate(a, color=color) for a, color in zip(arrows, sticker_colors(CORNER))))

        edge_arrows = sticker_arrows(EDGE)
        edge_matrix = matrix_tex(np.diag(EDGE), column_colors(EDGE)).move_to(matrix)
        edge_label = self.hud(address_tex(EDGE, font_size=40).to_corner(UL))
        with self.voice("Now try an edge: c equals one, zero, one. Its middle coordinate is zero, so "
                        "the middle column is all zeros. This cubelet has no sticker facing along "
                        "y, and in the matrix, that sticker simply drops out."):
            self.play(FadeOut(arrows), FadeOut(names), FadeOut(decomposition), FadeIn(edge_label),
                      cube.animate.ghosted([EDGE], 0.1), run_time=1)
            self.hud(edge_matrix)
            self.play(Draw(edge_arrows), ReplacementTransform(matrix, edge_matrix), run_time=1.5)
            self.when_said("middle column")
            self.play(Indicate(edge_matrix.get_columns()[1], color=GREY_B, scale_factor=1.4))

        center_arrows = sticker_arrows(CENTER)
        center_matrix = matrix_tex(np.diag(CENTER), column_colors(CENTER)).move_to(edge_matrix)
        center_label = self.hud(address_tex(CENTER, font_size=40).to_corner(UL))
        rank = self.hud(VGroup(
            MathTex(r"\#\,\text{non-zero columns} = \mathrm{rank}\,\mathrm{diag}(c)", font_size=32),
            MathTex(r"= \|c\|_1 = \#\,\text{stickers}", font_size=32),
        ).set_color(ACCENT).arrange(DOWN, aligned_edge=RIGHT).next_to(diag_group, DOWN, 0.5)
            .to_edge(RIGHT))  # Clear of the cube's right face.
        with self.voice("A center has just one non-zero column. So the number of non-zero "
                        "columns, the rank of the matrix, is once again the sticker count."):
            self.hud(center_matrix)
            # The center's one sticker arrow runs along the z axis, which steps aside for it.
            self.play(FadeOut(edge_arrows), cube.animate.ghosted([CENTER], 0.1),
                      FadeOut(axes[0][2]), FadeOut(axes[1][2]),
                      ReplacementTransform(edge_label, center_label),
                      ReplacementTransform(edge_matrix, center_matrix), run_time=1.2)
            self.play(Draw(center_arrows))
            self.when_said("So the number")
            self.play(Write(rank), run_time=2)

        bullets = self.hud(VGroup(*(MathTex(r"\bullet\;\;", r"\text{%s}" % head, r"\text{: %s}" % body,
                                            font_size=36)
                                    for head, body in [("how many stickers", "its rank"),
                                                       ("which way each faces", "its non-zero columns"),
                                                       ("what color each is", "the very same columns")]))
                           .arrange(DOWN, aligned_edge=LEFT, buff=0.4).move_to(1.2 * DOWN))
        for bullet in bullets:
            bullet[1].set_color(ACCENT)
        with self.voice("So this one little matrix tells us three things about a cubelet's "
                        "stickers."):
            self.play(*(FadeOut(m) for m in self.mobjects if m is not center_matrix and
                        m is not diag), run_time=1.2)
            self.play(VGroup(diag, center_matrix).animate.move_to(1.6 * UP))
        for bullet, words in zip(bullets, ["How many there are: its rank.",
                                           "Which way each one faces: its non-zero columns.",
                                           "And what color each one is: those very same columns."]):
            with self.voice(words, pause=0.6):
                self.play(FadeIn(bullet, shift=0.3 * RIGHT), run_time=0.8)
        self.wait(1)


# --- 7. A cubelet's configuration (c, R) ----------------------------------------

def twisted_corner_tex():
    """The twisted corner, in symbols: its position is home, its stickers are not."""
    twisted = VGroup(MathTex(r"R\,c = c", r"\;\checkmark", font_size=36),
                     MathTex(r"R\,\mathrm{diag}(c) \neq \mathrm{diag}(c)", r"\;\times",
                             font_size=36)).arrange(DOWN, aligned_edge=LEFT)
    twisted[0][1].set_color(GREEN)
    twisted[1][1].set_color(RED)
    return twisted


def product_tex(left, right, result, left_colors, right_colors, result_colors):
    """left · right = result, as bracketed integer matrices."""
    return VGroup(matrix_tex(left, left_colors), matrix_tex(right, right_colors),
                  MathTex("="), matrix_tex(result, result_colors)).arrange(RIGHT, buff=0.15)


class Configuration(Narrated):
    def construct(self):
        camera = dict(phi=64 * DEGREES, theta=22 * DEGREES)
        self.set_camera_orientation(**camera, zoom=0.92)
        axes = self.show_axes()
        cube = CubeMobject()
        self.add(cube)
        corner = cube.pieces[CORNER]
        with self.voice("Now let's scramble. As moves pile up, a cubelet gets carried around the "
                        "cube. But it's one solid little block: all it can do is move and turn "
                        "as a whole. And every move turns it about the origin."):
            self.play(cube.animate.ghosted([CORNER], 0.12), run_time=1.2)
            arrows = sticker_arrows(CORNER)
            self.play(Draw(arrows))
            corner.add(arrows)
            self.when_said("moves pile up")
            # Every move carries our corner somewhere new, and none brings it home: home and
            # twisted is the reveal for later in this chapter.
            scramble = [TURN_TOP, ((0, -1, 0), 1), (FRONT, -1)]
            for move in scramble:
                self.play(cube.turn(move), run_time=1.3)
                assert eigencube.position(CORNER, dict(cube.state)[CORNER]) != CORNER

        rule = self.hud(rotation_rule(font_size=28))
        rule_box = self.hud(SurroundingRectangle(rule, color=ACCENT, buff=0.15), overlay=True)
        VGroup(rule, rule_box).to_corner(UR)
        with self.voice("Remember our rule from the review: a rotation about the origin is a "
                        "matrix. And one turn after another is still just a rotation. So however "
                        "a cubelet got to where it is, its current orientation is a single "
                        "rotation matrix, which we call R."):
            self.when_said("a rotation about")
            self.play(Draw(rule_box), Write(rule), run_time=2)
        with self.voice("Let's put our corner back where it started, and look at it more closely."):
            self.undo(cube, scramble, run_time=0.8)

        identity = np.identity(3, dtype=int)
        config_title = self.hud(MathTex(r"\text{configuration: }(c,\, R)", font_size=40,
                                        color=ACCENT).to_corner(UL))
        c_tex = self.hud(address_tex(CORNER, font_size=36, color=POSITION)
                         .next_to(config_title, DOWN, 0.3, aligned_edge=LEFT))
        r_label = self.hud(labeled_matrix("R =", identity, BASIS_COLORS).scale(0.85)
                           .next_to(c_tex, DOWN, 0.3, aligned_edge=LEFT))
        with self.voice("That gives us the full description of a cubelet: its configuration, the "
                        "pair c and R. c names the cubelet by its home address. The matrix R "
                        "tells us how it's been turned. In the solved cube, every cubelet's R is the identity "
                        "matrix: no turn at all."):
            self.when_said("its configuration")
            self.play(Write(config_title))
            self.when_said("c names")
            self.play(FadeIn(c_tex))
            self.when_said("The matrix R tells")
            self.play(FadeIn(r_label))

        quarter = eigencube.rotation_matrix(TURN_TOP)
        r_top = self.hud(labeled_matrix("R =", quarter, BASIS_COLORS).scale(0.85)
                         .move_to(r_label, aligned_edge=LEFT))
        with self.voice("Let's turn the top face. For our corner, R becomes this quarter turn."):
            self.play(cube.turn(TURN_TOP), run_time=2)
            self.when_said("R becomes")
            self.play(ReplacementTransform(r_label, r_top))

        home = position_arrow(CORNER)
        moved = position_arrow(CORNER)  # Turns into p below, the way R turns c.
        p_eq = self.hud(VGroup(MathTex(r"p = R\,c =", font_size=36),
                               product_tex(quarter, [[x] for x in CORNER],
                                           [[x] for x in quarter @ CORNER],
                                           BASIS_COLORS, [POSITION], [POSITION]).scale(0.7))
                        .arrange(RIGHT).to_corner(DR))  # Clear of the x axis, lower left.
        above_captions(p_eq)
        home_label = self.facing_camera(MathTex("c", color=POSITION, font_size=40)
                                        .add_background_rectangle(opacity=0.75)
                                        .move_to(at(CORNER) + 0.55 * OUT + 0.3 * LEFT))
        moved_label = self.facing_camera(MathTex(r"p = R\,c", color=POSITION, font_size=40)
                                         .add_background_rectangle(opacity=0.75)
                                         .move_to(at((1, -1, 1)) + 0.8 * OUT + 1.15 * DOWN))
        with self.voice("So where is our cubelet, c, now? Take its home address, and rotate it: "
                        "apply the matrix R to c. We call the result p, the cubelet's position.",
                        pause=0.8):
            # From the front, c and p are equally foreshortened; from the usual angle, c points
            # almost straight at the camera and would show as a stub.
            # The x axis, seen end-on from there, steps aside meanwhile.
            self.play(corner.animate.set_opacity(0.45),  # See-through, so the arrow shows.
                      FadeOut(axes[0][0]), FadeOut(axes[1][0]))
            self.move_camera(phi=74 * DEGREES, theta=0, run_time=1.5)
            self.play(Draw(home), FadeIn(home_label))
            self.when_said("rotate it")
            self.add(moved)
            self.play(home.animate.set_opacity(0.35),
                      Rotate(moved, about_point=ORIGIN, **move_angle_axis(TURN_TOP)), run_time=2)
            self.play(FadeIn(moved_label), Write(p_eq), run_time=1.5)
        with self.voice("One, minus one, one. Front, left, top: exactly where the corner went."):
            self.play(Indicate(p_eq[1][-1], color=POSITION), Indicate(moved, color=POSITION))

        sticker_directions = quarter @ np.diag(CORNER)
        s_eq = self.hud(VGroup(MathTex(r"R\,\mathrm{diag}(c) =", font_size=36),
                               matrix_tex(sticker_directions, BASIS_COLORS).scale(0.7))
                        .arrange(RIGHT).move_to(p_eq, aligned_edge=LEFT))
        with self.voice("And where do its stickers point? Same idea: apply the matrix R to diag of c. "
                        "Matrix multiplication works column by column, so this rotates every "
                        "sticker at once. Green now points left, and red points to the front."):
            self.play(FadeOut(home), FadeOut(moved), FadeOut(home_label), FadeOut(moved_label),
                      FadeOut(p_eq), corner.animate.set_opacity(1))
            self.move_camera(**camera, run_time=1.5)
            self.play(FadeIn(axes[0][0]), FadeIn(axes[1][0]))
            self.when_said("apply the matrix R to diag")
            self.play(Write(s_eq), run_time=2)
            self.when_said("Green now")
            self.play(Indicate(s_eq[1].get_columns()[0], color=X_COLOR))
            self.when_said("and red")
            self.play(Indicate(s_eq[1].get_columns()[1], color=Y_COLOR))

        r_twist = self.hud(labeled_matrix("R =", dict(eigencube.apply_move_to_cube(
            TURN_FRONT, cube.state))[CORNER], BASIS_COLORS).scale(0.85)
            .move_to(r_top, aligned_edge=LEFT))
        twisted = self.hud(above_captions(twisted_corner_tex().to_corner(DR)))
        with self.voice("Turn the front face, and R picks up another factor. The corner is back "
                        "home: R c equals c. But look: its stickers are twisted. So the matrix R, "
                        "applied to diag of c, does not give back diag of c."):
            self.play(FadeOut(s_eq))
            self.play(cube.turn(TURN_FRONT), run_time=2)
            assert_twisted_home(cube.state)
            self.play(ReplacementTransform(r_top, r_twist))
            self.when_said("R c equals c")
            self.play(FadeIn(twisted[0]))
            self.when_said("But look")
            self.play(*(Indicate(a, color=color) for a, color in zip(arrows, sticker_colors(CORNER))))
            self.when_said("does not give back")
            self.play(FadeIn(twisted[1]))

        code = self.hud(above_captions(code_listing("solved_cube").scale(0.62).to_corner(DL)))
        with self.voice("And that is the complete state of the cube: one configuration per "
                        "cubelet. No list of stickers, no flags for twisted corners or flipped "
                        "edges. Nothing else."):
            self.play(FadeOut(twisted))
            self.when_said("one configuration")
            self.play(FadeIn(code, shift=0.2 * UP))
        self.wait(1)


# --- 8. Moves: a dot product and a matrix product -------------------------------

class Moves(Narrated):
    def construct(self):
        theta = 26 * DEGREES
        self.set_camera_orientation(phi=66 * DEGREES, theta=theta, zoom=1.1)
        axes = self.show_axes()
        cube = CubeMobject(eigencube.shuffle(eigencube.solved_cube, 12, seed=8))
        self.add(cube)
        self.say("So how do we make a move? Say we turn the top face. Which cubelets does that "
                 "turn?")

        v_arrow = arm(TOP, 2.6, ACCENT, thickness=0.05)
        v_label = self.hud(MathTex(r"\mathbf{v} = \mathbf{e}_z", color=ACCENT, font_size=48)
                           .to_corner(UR))
        with self.voice("The ones in the top layer, with z equal to one. To say that for any face, "
                        "take the face's axis, v. For the top face, v is e z, pointing up."):
            self.play(cube.animate.ghosted(cube.slice_cubelets(TURN_TOP), 0.15), run_time=1.5)
            self.when_said("the face's axis")
            # The z axis steps aside, so that v, which runs along it, stands out.
            self.play(FadeOut(axes[0][2]), FadeOut(axes[1][2]), Draw(v_arrow), FadeIn(v_label))

        dot = self.hud(VGroup(
            MathTex(r"\mathbf{v} \cdot p", r"= \text{how far } p \text{ reaches along } \mathbf{v}",
                    font_size=36),
            Text("Essence of Linear Algebra,\nch. 9: Dot products and duality", font_size=20,
                 color=ACCENT, line_spacing=0.8),
        ).arrange(DOWN, aligned_edge=LEFT).to_corner(UL))
        derivation = self.hud(MathTex(r"\mathbf{v} \cdot p &= %s \cdot (x,\, y,\, z) \\ "
                                      r"&= 0x + 0y + 1z = z" % vector_tex_string(TOP), font_size=36)
                              .next_to(v_label, DOWN, 0.4, aligned_edge=RIGHT))
        # Each label sits level with its layer: at the middle of the layer's stripe on the right
        # face (perspective spaces the face's front edge wider than its back edge), moved along the
        # camera's horizontal to just past the cube.
        right = screen_right(theta)
        to_past_cube = at((-1.5, 0, 0)) @ right + 0.55  # From the stripe's middle, sideways.
        layer_labels = [self.facing_camera(MathTex(text, color=ACCENT, font_size=40)
                                           .move_to(at((0, 1.5, z)) + to_past_cube * right))
                        for z, text in [(1, "+1"), (0, "0"), (-1, "-1")]]
        with self.voice("Now take the dot product of v with each cubelet's current position, p. "
                        "Dotting with a unit axis measures how far p reaches along it. Essence of "
                        "Linear Algebra has a whole chapter on why that works: dot products and "
                        "duality."):
            self.play(cube.animate.ghosted(cube.pieces), run_time=1)
            self.when_said("Dotting")
            self.play(Write(dot[0]), run_time=2)
            self.when_said("Essence of")
            self.play(FadeIn(dot[1]))
        with self.voice("To compute it, multiply matching coordinates and add up: zero times x, "
                        "plus zero times y, plus one times z. That's simply z, the height."):
            self.play(Write(derivation), run_time=3)
        with self.voice("That's one for the top layer, zero for the middle, and minus one for the "
                        "bottom."):
            self.play(FadeOut(axes[0][1]), FadeOut(axes[1][1]))  # The y axis runs where "0" goes.
            for label, words in zip(layer_labels, ["one for", "zero for", "minus one"]):
                self.when_said(words)
                self.play(FadeIn(label), run_time=0.6)

        selector = self.hud(MathTex(r"\mathbf{v} \cdot (R\,c) > 0", font_size=48)
                            .next_to(dot, DOWN, 0.4, aligned_edge=LEFT))
        # Tiled, because Cairo depth-sorts each shape as a whole by its center.
        plane = VGroup(*(Square(0.9 * SPACING, fill_color=ACCENT, fill_opacity=0.2, stroke_width=0)
                         for _ in range(16))).arrange_in_grid(4, 4, buff=0)
        plane.set_shade_in_3d(True).shift(0.5 * SPACING * OUT)
        with self.voice("So the top layer is exactly where v dot p is positive. And since p is R c, "
                        "the test is: v dot R c, greater than zero. One dot product per cubelet, "
                        "and the layer selects itself."):
            self.play(cube.animate.ghosted(cube.slice_cubelets(TURN_TOP), 0.15), FadeIn(plane))
            self.when_said("the test is")
            self.play(Write(selector), run_time=2)
        self.play(FadeOut(plane), FadeOut(*layer_labels), FadeOut(dot), FadeOut(derivation),
                  selector.animate.to_corner(UL))
        self.play(FadeIn(axes[0][1]), FadeIn(axes[1][1]), run_time=0.6)  # Once "0" is gone.

        turn = eigencube.rotation_matrix(TURN_TOP)
        basis = VGroup(*(overlay_arrow(ORIGIN, at(unit(i), 2.0), color, thickness=0.05)
                         for i, color in enumerate(BASIS_COLORS)))
        m_group = self.hud(labeled_matrix("M =", turn, BASIS_COLORS).to_corner(UR))
        with self.voice("Next, the turn itself: a quarter rotation, M. To write it down, use the "
                        "trick from the review. Watch where the basis vectors land: those are the "
                        "columns."):
            # The axes step aside, so that the turning basis vectors are not confused with them.
            # The basis is drawn over everything, so the cube turns see-through first: drawn
            # meanwhile, the arrows would cover the top layer's still-opaque stickers.
            # Only the x and y axes: the z axis already stepped aside for v, and fading it out
            # again would first put it back, over the opaque top layer.
            self.play(FadeOut(v_arrow), FadeOut(v_label), *(FadeOut(axes[part][i])
                                                            for part in (0, 1) for i in (0, 1)),
                      cube.animate.ghosted(opacity=0.07), run_time=0.6)
            self.play(Draw(basis), run_time=0.8)
            self.play(Write(m_group[0]), FadeIn(m_group[1].get_brackets()))
            self.when_said("Watch")
            self.play(Rotate(basis, **move_angle_axis(TURN_TOP), about_point=ORIGIN), run_time=4)
        # Each label just past its arrow's tip: above the flat ones, beside the upright one.
        landings = [self.facing_camera(vector_tex(turn[:, i], color=color, font_size=34)
                                       .add_background_rectangle(opacity=0.75)
                                       .move_to(at(turn[:, i], 2.6) + 0.75 * OUT if turn[2, i] == 0
                                                else at(turn[:, i], 2.0) + 0.75 * right))
                    for i, color in enumerate(BASIS_COLORS)]
        for column, a, color, landing, (words, cue) in zip(
                m_group[1].get_columns(), basis, BASIS_COLORS, landings,
                [("e x lands on minus e y: zero, minus one, zero.", "lands on"),
                 ("e y lands on e x: one, zero, zero.", "lands on"),
                 ("And e z stays put: zero, zero, one. It's the axis of the turn: our eigenvector "
                  "again.", "stays put")]):
            with self.voice(words, pause=0.25):
                self.when_said(cue)
                self.play(FadeIn(landing), Indicate(a, color=color), run_time=0.8)
                self.play(FadeIn(column, shift=0.3 * DOWN), run_time=0.8)

        update = self.hud(MathTex(r"R \;\leftarrow\; M\,R", font_size=52, color=ACCENT).next_to(
            selector, DOWN, 0.5, aligned_edge=LEFT))
        with self.voice("Then every selected cubelet multiplies its R by M: R becomes M times R. "
                        "That's the whole move."):
            self.play(FadeOut(basis), FadeOut(*landings), FadeIn(axes[0]), FadeIn(*axes[1]),
                      cube.animate.ghosted(cube.pieces))
            self.when_said("R becomes")
            self.play(Write(update))
            self.play(cube.turn(TURN_TOP), run_time=2)

        with self.voice("One dot product to select, one matrix product to turn. And since every "
                        "entry is zero, one, or minus one, the arithmetic stays exact, no matter "
                        "how many moves you make."):
            self.play(FadeOut(m_group))  # It's the top turn's matrix; other faces turn next.
            self.turn_while_speaking(cube, [TURN_FRONT, ((0, 1, 0), 1), ((0, 0, -1), -1),
                                            ((-1, 0, 0), 1)], run_time=0.9)
        code = self.hud(above_captions(code_listing("apply_move_to_cubelet_rotation").scale(0.65)
                                       .to_corner(DL)))
        with self.voice("Here it is in the actual code. That really is the entire move logic."):
            # The cube fades back, so that the wide listing can overlap it and stay legible.
            self.play(FadeIn(code, shift=0.2 * UP), cube.animate.ghosted(opacity=0.25),
                      FadeOut(axes[0]), FadeOut(*axes[1]))
        self.wait(1)


# --- 9. When is it solved? ------------------------------------------------------

class Solved(Narrated):
    def construct(self):
        theta = 22 * DEGREES
        self.set_camera_orientation(phi=64 * DEGREES, theta=theta, zoom=1.05)
        axes = self.show_axes()
        state = eigencube.apply_move_to_cube(TURN_FRONT, eigencube.apply_move_to_cube(
            TURN_TOP, eigencube.solved_cube))
        cube = CubeMobject(state)
        self.add(cube)
        assert_twisted_home(state)
        rotation = dict(state)[CORNER]
        self.play(cube.animate.ghosted([CORNER], 0.12))
        arrows = sticker_arrows(CORNER).apply_matrix(np.array(rotation, dtype=float),
                                                      about_point=ORIGIN)
        with self.voice("Last question: when is the cube solved? A cubelet is solved when each of "
                        "its stickers faces the same way as the center of its color."):
            self.play(Draw(arrows))
        criterion = self.hud(MathTex(r"R\,\mathrm{diag}(c)", r"=", r"\mathrm{diag}(c)",
                                     font_size=52).to_corner(UL))
        with self.voice("In symbols: apply the matrix R to diag of c, and you should get diag of "
                        "c right back."):
            self.when_said("apply")
            self.play(Write(criterion), run_time=2)

        twisted = self.hud(twisted_corner_tex().next_to(criterion, DOWN, 0.5, aligned_edge=LEFT))
        with self.voice("Remember our twisted corner? Its position is right: R c equals c. But its "
                        "stickers are not, and the test catches it."):
            self.when_said("Its position")
            self.play(FadeIn(twisted[0]))
            self.when_said("But its")
            self.play(FadeIn(twisted[1]))
        self.play(FadeOut(twisted), FadeOut(arrows))

        with self.voice("Here's a lovely subtlety. For a corner, three stickers facing home force R "
                        "to be the identity. The same goes for an edge: pin down two directions, "
                        "and a rotation has nowhere left to go."):
            self.play(cube.turn((FRONT, -1)), run_time=1.2)
            self.play(cube.turn((TOP, -1)), run_time=1.2)
            self.when_said("force R")
            self.play(Indicate(criterion, color=ACCENT))

        with self.voice("But a center has just one sticker. We drew a little mark on this one, so "
                        "that you can watch it turn."):
            # Clear view of the mark: no axes on top of it, and a camera looking down on it.
            self.play(cube.animate.ghosted([CENTER], 0.12), FadeOut(axes[0]), FadeOut(*axes[1]))
            # The cube also steps right, to make room for the calculation that follows.
            self.move_camera(phi=30 * DEGREES, theta=theta, zoom=1.05,
                             frame_center=-3.4 * screen_right(theta), run_time=1.5)
            self.when_said("watch it turn")
            self.play(cube.turn(TURN_TOP), run_time=2.5)
        spun = np.array(dict(cube.state)[CENTER])
        r_center = labeled_matrix("R =", spun, BASIS_COLORS).scale(0.8)
        product = VGroup(MathTex(r"R\,\mathrm{diag}(c) =", font_size=32),
                         product_tex(spun, np.diag(CENTER), spun @ np.diag(CENTER), BASIS_COLORS,
                                     column_colors(CENTER), column_colors(CENTER)).scale(0.6),
                         MathTex(r"\checkmark", font_size=32))
        product[2].set_color(GREEN)
        worked = self.hud(VGroup(r_center, product.arrange(RIGHT))
                          .arrange(DOWN, aligned_edge=LEFT, buff=0.3)
                          .next_to(criterion, DOWN, 0.5, aligned_edge=LEFT))
        with self.voice("Its rotation is a quarter turn, not the identity. But apply it to diag of "
                        "c: only the last column is non-zero, and the turn leaves it in place. So "
                        "we get diag of c right back, and the test passes."):
            self.play(FadeIn(r_center))
            self.when_said("But apply")
            self.play(Write(product[0]), FadeIn(product[1]), run_time=1.5)
            self.when_said("we get diag")
            self.play(Write(product[2]))
        self.say("On a plain cube, that spin is invisible, and the encoding, rightly, doesn't care "
                 "either. It's exactly as picky as the colors are.")

        code = self.hud(above_captions(code_listing("is_cubelet_solved").scale(0.7)
                                       .to_corner(DL)))
        with self.voice("And in code, the whole test is three lines, straight from the math."):
            # The worked example makes room; the criterion stays, next to its code.
            self.play(FadeOut(worked), FadeIn(code, shift=0.2 * UP))
        self.wait(1)


# --- 10. Recap and outro ----------------------------------------------------------

class Outro(Narrated):
    def construct(self):
        self.wait(1.5)  # A breath before the recap (which the README's summary table mirrors).
        recap = VGroup(*(MathTex(r"\text{%s}" % name, r"\;%s\;" % relation, rhs, font_size=36)
                         for name, relation, rhs in [
            ("cubelets", "=", r"\text{vectors } c \in \{-1, 0, 1\}^3"),
            ("cubelet type", "=", r"\|c\|_1 = \#\,\text{stickers}"),
            ("colors", "=", r"\pm\mathbf{e}_x,\ \pm\mathbf{e}_y,\ \pm\mathbf{e}_z"),
            ("sticker colors", "=", r"\text{columns of } \mathrm{diag}(c)"),
            ("state", "=", r"\text{one configuration } (c,\, R) \text{ per cubelet}"),
            ("position", "=", r"p = R\,c"),
            ("sticker directions", "=", r"R\,\mathrm{diag}(c)"),
            ("move", ":", r"\mathbf{v}\cdot(R\,c) > 0 \;\Rightarrow\; R \leftarrow M R"),
            ("solved", r"\iff", r"R\,\mathrm{diag}(c) = \mathrm{diag}(c)"),
        ])).arrange(DOWN, aligned_edge=LEFT, buff=0.2).move_to(0.35 * UP)
        self.label(recap)
        for row in recap:
            row[0].set_color(ACCENT)
        self.say("Let's step back, and look at what we built.")
        for row, words in zip(recap, [
                "Cubelets are vectors, c, with coordinates minus one, zero, or one.",
                "The Manhattan norm of c tells us its type: how many stickers it has.",
                "Colors are basis vectors.",
                "A cubelet's sticker colors are the columns of diag of c.",
                "The state of the cube is one configuration, c and R, per cubelet.",
                "Where a cubelet is now: its position p, the matrix R applied to c.",
                "Which way its stickers face now: the matrix R applied to diag of c.",
                "A move is a dot product to select, and a matrix product to turn.",
                "And a cubelet is solved when the matrix R, applied to diag of c, changes nothing."]):
            with self.voice(words, pause=0.4):
                self.play(FadeIn(row, shift=0.3 * RIGHT), run_time=0.9)
        with self.voice("That's the entire model. Everything else falls out of the geometry."):
            self.play(Circumscribe(recap, color=ACCENT, buff=0.25), run_time=1.5)
        self.play(FadeOut(recap))

        self.set_camera_orientation(**{**CAMERA, "zoom": 1.3})
        scrambled = eigencube.shuffle(eigencube.solved_cube, 40, seed=11)
        cube = CubeMobject(scrambled)
        solution = eigencube.solve(scrambled)
        with self.voice("On top of that model, Eigencube solves the cube layer by layer: an A-star "
                        "search, guided by how many quarter turns each cubelet is from home, plus "
                        "one classic move sequence to finish."):
            self.play(FadeIn(cube, scale=0.8))
            self.begin_ambient_camera_rotation(rate=0.15)
        with self.voice("Here it is, solving a scramble."):
            for move in solution:
                self.add_sound(str(click()), gain=-14)  # Otherwise the solve plays in silence.
                self.play(cube.turn(move), run_time=0.09, rate_func=linear)
        assert eigencube.is_cube_solved(cube.state)
        self.stop_ambient_camera_rotation()
        self.move_camera(**CAMERA, run_time=1.5)
        self.wait(1.5)  # Let the solved cube land.

        logo = SVGMobject(str(HERE.parent / "img" / "logo.svg")).scale(0.8)
        card = self.hud(VGroup(
            logo,
            Text("Watch: Essence of Linear Algebra", font_size=40, weight=BOLD),
            VGroup(Text("3Blue1Brown", font_size=26, color=ACCENT),
                   Text("youtube.com/@3blue1brown", font="DejaVu Sans Mono", font_size=22,
                        color=ACCENT)).arrange(RIGHT, buff=0.4),
            Text("Code: github.com/smolkaj/eigencube", font="DejaVu Sans Mono", font_size=24,
                 color=GREY_B),
        ).arrange(DOWN, buff=0.3).move_to(0.5 * UP))
        music_start = self.renderer.time
        with self.voice("If this way of seeing made you smile, the credit goes to Essence of Linear "
                        "Algebra. Go watch it. You'll start seeing matrices as motions everywhere, "
                        "even inside a puzzle from 1974.", pause=1.4):
            self.play(FadeOut(cube), run_time=1)
            self.play(FadeIn(card, shift=0.3 * UP), run_time=1.5)
        # The last word calls back to the opening's sticker list.
        self.say("It was never a list of fifty-four colors. It was linear algebra, all along.")
        self.wait(12)  # The music plays on under the end card before it rings out.
        self.add_music_since(music_start)


SCENES = [StickerNightmare, LinearAlgebraReview, FixedFrame, CountingStickers, ColorsAreVectors,
          DiagTrick, Configuration, Moves, Solved, Outro]
