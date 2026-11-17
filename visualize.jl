using GLMakie
using GeometryBasics: Mesh, Point3f, TriangleFace
include(joinpath(@__DIR__, "src", "TorusGeometry.jl"))
using .TorusGeometry

const TORUS = Torus(2.0, 0.7)
const SEED = seed_angles(42)
const SEGMENTS = 48

point3(p) = Point3f(p[1], p[2], p[3])

function line_points(edges; loop = false)
	result = Point3f[]
	for edge in edges
		append!(result, point3.(edge))
		loop && push!(result, point3(first(edge)))
		push!(result, Point3f(NaN, NaN, NaN))
	end
	return result
end

function frame(time)
	state = growth(time)
	patch = integrate_patch(TORUS, SEED, state.u_span, state.v_span, SEGMENTS)
	data = patch_mesh(patch; u_closed = state.u_closed, v_closed = state.v_closed)
	geometry = Mesh(point3.(data.vertices), [TriangleFace(face...) for face in data.faces])
	colors = [area_scale(TORUS, u, v) for (u, v) in data.parameters]
	u_edges = line_points(boundary_points(data, :u, state.u_closed); loop = state.v_closed)
	v_edges = line_points(boundary_points(data, :v, state.v_closed); loop = state.u_closed)
	return (; state, geometry, colors, u_edges, v_edges)
end

initial = frame(0.0)
geometry = Observable(initial.geometry)
colors = Observable(initial.colors)
u_edges = Observable(initial.u_edges)
v_edges = Observable(initial.v_edges)
parameter_us = Observable(collect(range(SEED[1] - initial.state.u_span / 2,
									  SEED[1] + initial.state.u_span / 2; length = SEGMENTS + 1)))
parameter_vs = Observable(collect(range(SEED[2] - initial.state.v_span / 2,
									  SEED[2] + initial.state.v_span / 2; length = SEGMENTS + 1)))
parameter_area_scale = Observable([area_scale(TORUS, u, v) for u in parameter_us[], v in parameter_vs[]])
parameter_u_edges = Observable(Point2f[])
parameter_v_edges = Observable(Point2f[])
parameter_grid = Observable(Point2f[])
caption = Observable("")

fig = Figure(size = (1400, 800))
Label(fig[1, 1:3], "Jacobian Torus: integrate local tangents, identify boundaries", fontsize = 22)
axis = Axis3(fig[2, 1], aspect = :data, xlabel = "x", ylabel = "y", zlabel = "z")
limits!(axis, -3.0, 3.0, -3.0, 3.0, -1.0, 1.0)
plot = mesh!(axis, geometry; color = colors, colormap = :thermal,
			 colorrange = (TORUS.r * (TORUS.R - TORUS.r), TORUS.r * (TORUS.R + TORUS.r)),
			 shading = false)
Colorbar(fig[2, 3], plot, label = "Local area scale |x_u × x_v|")
lines!(axis, u_edges; color = :orange, linewidth = 4)
lines!(axis, v_edges; color = :cyan, linewidth = 4)
parameter_axis = Axis(fig[2, 2], title = "Parameter domain: Jacobian area scale", xlabel = "u", ylabel = "v", aspect = DataAspect())
limits!(parameter_axis, SEED[1] - pi - 0.3, SEED[1] + pi + 0.3,
		SEED[2] - pi - 0.3, SEED[2] + pi + 0.3)
heatmap!(parameter_axis, parameter_us, parameter_vs, parameter_area_scale;
		 colormap = :thermal,
		 colorrange = (TORUS.r * (TORUS.R - TORUS.r), TORUS.r * (TORUS.R + TORUS.r)))
lines!(parameter_axis, parameter_grid; color = :gray, linewidth = 1)
lines!(parameter_axis, parameter_u_edges; color = :orange, linewidth = 4)
lines!(parameter_axis, parameter_v_edges; color = :cyan, linewidth = 4)
scatter!(parameter_axis, [Point2f(SEED...)]; color = :black, markersize = 10)
seed_point = torus_point(TORUS, SEED...)
xu, xv = tangents(TORUS, SEED...)
scatter!(axis, [point3(seed_point)]; color = :white, markersize = 12)
lines!(axis, point3.([seed_point, seed_point + 0.3xu]); color = :red, linewidth = 5)
lines!(axis, point3.([seed_point, seed_point + 0.3xv]); color = :blue, linewidth = 5)
Label(fig[3, 1:3], caption, tellwidth = false)
controls = fig[4, 1:3] = GridLayout()
slider = Slider(controls[1, 1], range = 0.0:0.005:1.0, startvalue = 0.0)
button = Button(controls[1, 2], label = "Play / pause")
playing = Ref(false)
play_generation = Ref(0)

function update(time)
	current = frame(time)
	colors[] = current.colors
	geometry[] = current.geometry
	u_edges[] = current.u_edges
	v_edges[] = current.v_edges
	umin, umax = SEED[1] .+ (-current.state.u_span / 2, current.state.u_span / 2)
	vmin, vmax = SEED[2] .+ (-current.state.v_span / 2, current.state.v_span / 2)
	parameter_us[] = collect(range(umin, umax; length = SEGMENTS + 1))
	parameter_vs[] = collect(range(vmin, vmax; length = SEGMENTS + 1))
	parameter_area_scale[] = [area_scale(TORUS, u, v) for u in parameter_us[], v in parameter_vs[]]
	gap = Point2f(NaN, NaN)
	parameter_u_edges[] = [Point2f(umin, vmin), Point2f(umin, vmax), gap,
						   Point2f(umax, vmin), Point2f(umax, vmax)]
	parameter_v_edges[] = [Point2f(umin, vmin), Point2f(umax, vmin), gap,
						   Point2f(umin, vmax), Point2f(umax, vmax)]
	grid = Point2f[]
	for fraction in (0.25, 0.5, 0.75)
		u = umin + fraction * (umax - umin)
		v = vmin + fraction * (vmax - vmin)
		append!(grid, [Point2f(u, vmin), Point2f(u, vmax), gap,
					   Point2f(umin, v), Point2f(umax, v), gap])
	end
	parameter_grid[] = grid
	u_status = current.state.u_closed ? "identified" : "open"
	v_status = current.state.v_closed ? "identified" : "open"
	caption[] = "t = $(round(time; digits = 3))   |   orange u: $u_status   |   cyan v: $v_status   |   red x_u, blue x_v"
end
on(slider.value) do time
	update(time)
end
update(0.0)

on(button.clicks) do _
	playing[] = !playing[]
	play_generation[] += 1
	generation = play_generation[]
	if playing[]
		slider.value[] >= 1 && set_close_to!(slider, 0.0)
		@async begin
			while playing[] && play_generation[] == generation && slider.value[] < 1
				sleep(0.05)
				(playing[] && play_generation[] == generation) || break
				set_close_to!(slider, min(slider.value[] + 0.005, 1.0))
			end
			play_generation[] == generation && (playing[] = false)
		end
	end
end

screen = display(fig)
wait(screen)
