using Test
using LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "TorusGeometry.jl"))
using .TorusGeometry

@testset "Jacobian torus" begin
	t = Torus(2.0, 0.7)
	@test_throws ArgumentError Torus(1.0, 1.0)
	@test seed_angles(42) == seed_angles(42)
	@test seed_angles(42) != seed_angles(43)
	seed = seed_angles(42)
	h = 1e-5
	for (u, v) in ((0.0, 0.0), (1.1, 2.2), (-2.0, 4.0))
		xu, xv = tangents(t, u, v)
		@test xu ≈ (torus_point(t, u + h, v) - torus_point(t, u - h, v)) / (2h) atol = 1e-9
		@test xv ≈ (torus_point(t, u, v + h) - torus_point(t, u, v - h)) / (2h) atol = 1e-9
		@test area_scale(t, u, v) ≈ t.r * (t.R + t.r * cos(v))
		@test cross(xu, xv) ⋅ [cos(u)*cos(v), sin(u)*cos(v), sin(v)] > 0
	end

	@test_throws ArgumentError growth(-0.1)
	@test_throws ArgumentError growth(1.1)
	@test growth(0).u_span > 0
	@test !growth(0).u_closed && !growth(0).v_closed
	@test growth(0.5).u_closed && !growth(0.5).v_closed
	@test growth(1).u_closed && growth(1).v_closed
	@test growth(1).u_span ≈ 2pi
	@test growth(1).v_span ≈ 2pi

	function max_error(n)
		patch = integrate_patch(t, seed, 2pi, 2pi, n)
		maximum(norm(patch.points[i, j] - torus_point(t, patch.us[i], patch.vs[j]))
				for i in eachindex(patch.us), j in eachindex(patch.vs))
	end
	@test max_error(48) < 1e-5
	@test max_error(24) < max_error(12) / 10
	@test_throws ArgumentError integrate_patch(t, seed, 1, 1, 3)

	patch = integrate_patch(t, seed, 2pi, 2pi, 8)
	for (uc, vc, vertices, faces, boundary_edges) in ((false, false, 81, 128, 32),
													   (true, false, 72, 128, 16),
													   (true, true, 64, 128, 0))
		mesh = patch_mesh(patch; u_closed = uc, v_closed = vc)
		@test length(mesh.vertices) == vertices
		@test length(mesh.faces) == faces
		edges = Dict{Tuple{Int,Int},Int}()
		for (a, b, c) in mesh.faces, (x, y) in ((a, b), (b, c), (c, a))
			edge = minmax(x, y)
			edges[edge] = get(edges, edge, 0) + 1
		end
		@test count(==(1), values(edges)) == boundary_edges
		@test all(n -> n == 1 || n == 2, values(edges))
		@test length(boundary_points(mesh, :u, uc)) == (uc ? 1 : 2)
		@test length(boundary_points(mesh, :v, vc)) == (vc ? 1 : 2)
	end
end
