module TorusGeometry

using LinearAlgebra
using Random

export Torus, seed_angles, torus_point, tangents, area_scale, growth,
	   integrate_patch, patch_mesh, boundary_points

struct Torus
	R::Float64
	r::Float64
	function Torus(R::Real = 2.0, r::Real = 0.7)
		0 < r < R || throw(ArgumentError("A ring torus requires R > r > 0"))
		new(R, r)
	end
end

seed_angles(seed::Integer = 42) = 2pi .* Tuple(rand(MersenneTwister(seed), 2))

torus_point(t::Torus, u, v) = [(t.R + t.r * cos(v)) * cos(u),
							   (t.R + t.r * cos(v)) * sin(u), t.r * sin(v)]

function tangents(t::Torus, u, v)
	rho = t.R + t.r * cos(v)
	xu = [-rho * sin(u), rho * cos(u), 0.0]
	xv = [-t.r * sin(v) * cos(u), -t.r * sin(v) * sin(u), t.r * cos(v)]
	return xu, xv
end

area_scale(t::Torus, u, v) = norm(cross(tangents(t, u, v)...))

function growth(time::Real; initial_span::Real = 0.18)
	0 <= time <= 1 || throw(ArgumentError("time must be between 0 and 1"))
	0 < initial_span < 2pi || throw(ArgumentError("initial span must be between 0 and 2pi"))
	smooth(x) = x * x * (3 - 2x)
	u = initial_span + (2pi - initial_span) * smooth(min(2time, 1))
	v = initial_span + (2pi - initial_span) * smooth(max(2time - 1, 0))
	return (u_span = u, v_span = v, u_closed = time >= 0.5, v_closed = time >= 1)
end

# Simpson quadrature integrates each local tangent over one mesh interval.
function step(t::Torus, fixed, a, b, direction::Int)
	derivative(s) = direction == 1 ? tangents(t, s, fixed)[1] : tangents(t, fixed, s)[2]
	return (b - a) / 6 .* (derivative(a) .+ 4 .* derivative((a + b) / 2) .+ derivative(b))
end

function integrate_patch(t::Torus, seed::Tuple{<:Real,<:Real}, u_span::Real, v_span::Real, segments::Integer)
	segments >= 2 && iseven(segments) || throw(ArgumentError("segments must be even and at least 2"))
	0 < u_span <= 2pi && 0 < v_span <= 2pi || throw(ArgumentError("spans must be in (0, 2pi]"))
	u0, v0 = seed
	us = collect(range(u0 - u_span / 2, u0 + u_span / 2; length = segments + 1))
	vs = collect(range(v0 - v_span / 2, v0 + v_span / 2; length = segments + 1))
	mid = segments ÷ 2 + 1
	points = Matrix{Vector{Float64}}(undef, segments + 1, segments + 1)
	points[mid, mid] = torus_point(t, u0, v0)
	for i in mid+1:segments+1
		points[i, mid] = points[i-1, mid] + step(t, v0, us[i-1], us[i], 1)
	end
	for i in mid-1:-1:1
		points[i, mid] = points[i+1, mid] + step(t, v0, us[i+1], us[i], 1)
	end
	for i in eachindex(us)
		for j in mid+1:segments+1
			points[i, j] = points[i, j-1] + step(t, us[i], vs[j-1], vs[j], 2)
		end
		for j in mid-1:-1:1
			points[i, j] = points[i, j+1] + step(t, us[i], vs[j+1], vs[j], 2)
		end
	end
	return (points = points, us = us, vs = vs)
end

function patch_mesh(patch; u_closed::Bool, v_closed::Bool)
	nu, nv = size(patch.points)
	u_closed && (nu -= 1)
	v_closed && (nv -= 1)
	vertices = [patch.points[i, j] for i in 1:nu for j in 1:nv]
	parameters = [(patch.us[i], patch.vs[j]) for i in 1:nu for j in 1:nv]
	index(i, j) = (i - 1) * nv + j
	faces = NTuple{3,Int}[]
	for i in 1:(u_closed ? nu : nu - 1), j in 1:(v_closed ? nv : nv - 1)
		next_i = i == nu ? 1 : i + 1
		next_j = j == nv ? 1 : j + 1
		a, b, c, d = index(i, j), index(next_i, j), index(next_i, next_j), index(i, next_j)
		push!(faces, (a, b, c), (a, c, d))
	end
	return (vertices = vertices, parameters = parameters, faces = faces, nu = nu, nv = nv)
end

function boundary_points(mesh, direction::Symbol, closed::Bool)
	index(i, j) = (i - 1) * mesh.nv + j
	if direction == :u
		first_edge = [mesh.vertices[index(1, j)] for j in 1:mesh.nv]
		second_edge = [mesh.vertices[index(mesh.nu, j)] for j in 1:mesh.nv]
	elseif direction == :v
		first_edge = [mesh.vertices[index(i, 1)] for i in 1:mesh.nu]
		second_edge = [mesh.vertices[index(i, mesh.nv)] for i in 1:mesh.nu]
	else
		throw(ArgumentError("direction must be :u or :v"))
	end
	return closed ? [first_edge] : [first_edge, second_edge]
end

end
julia --project=. .\test\runtests.jl