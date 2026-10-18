using TestItemRunner

@run_package_tests

@testitem "every move is a 4-cycle and its inverse undoes it" begin
    using EigencubePretty: moves, inverse
    c = scramble(1000, 7)
    for move in moves
        @test foldl((x, _) -> move(x), 1:4; init = c) == c
        @test move(c) != c
        @test inverse(move)(move(c)) == c
    end
end

@testitem "rotation matrices are proper rotations" begin
    using EigencubePretty: rotation_matrix, moves, identity3
    using LinearAlgebra: det
    for move in moves
        R = rotation_matrix(move)
        @test R' * R == identity3
        @test det(R) == 1
        @test R^4 == identity3
    end
end

@testitem "the random stream is the reference's xorshift64" begin
    using EigencubePretty: Xorshift64
    rng = Xorshift64(42)
    # x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17, starting from 42, then the top 53 bits.
    x = UInt64(42)
    x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
    @test rand(rng) === (x >> 11) * 0x1.0p-53
end

@testitem "solutions solve the cube" begin
    using EigencubePretty: apply_moves
    c = scramble(100_000, 1)
    @test !is_solved(c)
    # Solved means R·diag(c) = diag(c); centers may still be turned about their own axis.
    @test is_solved(apply_moves(c, solve(c)))
end

@testitem "solution lengths match the rust_pretty reference" begin
    # Solution lengths printed by eigencube's rust_pretty for seeds 1–3.
    @test [length(solve(scramble(100_000, s))) for s in 1:3] == [130, 150, 178]
end

@testitem "heuristics do not allocate" begin
    using EigencubePretty: top_layer_heuristic, middle_layer_heuristic, bottom_layer_edge_heuristic, bottom_layer_corner_heuristic
    allocated(f, c) = (f(c); @allocated f(c))
    c = scramble(1000, 3)
    for h in (top_layer_heuristic, middle_layer_heuristic, bottom_layer_edge_heuristic, bottom_layer_corner_heuristic)
        @test iszero(allocated(h, c))
    end
end

@testitem "an A* attempt allocates nothing once its buffers are sized" begin
    using EigencubePretty: Search, bounded_search!, heuristic, TopAndMiddleLayer, out_of_budget
    # An unreachable goal, so the attempt runs its whole budget and builds no path.
    s, c = Search(), scramble(100_000, 2)
    attempt() = bounded_search!(s, c, Returns(false), x -> heuristic(TopAndMiddleLayer(), 3, x), () -> 1.0, 20_000)
    @test attempt() === out_of_budget
    @test iszero(@allocated attempt())
end

@testitem "heuristic norm is (Σ √d)²" begin
    using EigencubePretty: pnorm_half
    @test pnorm_half(c -> 4, [1, 2], _ -> true) === 16.0
end
