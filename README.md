# Jacobian Torus: local geometry, global topology

This Julia prototype grows a small parameter-space patch into a torus. It is a learning experiment, separate from the empty Visual Studio C++ project; C++ and assembly can be explored later.

## Project thesis

From local Jacobian information to a realized global object. The Jacobian supplies local geometric information, including tangent directions and area scaling; integration realizes the surface, and boundary identification realizes its topology. Together, these structures realize the torus—and make the process visible through computation and animation.

## Run on this machine

Open PowerShell in `C:\Users\tmcla\JacobianTorus`. Julia 1.12 was found in the location below; substitute `julia` if Julia is on your PATH.

```powershell
$Julia = 'C:\Users\tmcla\AppData\Local\Programs\Julia-1.12.0\bin\julia.exe'
& $Julia --project=. -e 'using Pkg; Pkg.instantiate()'
& $Julia --project=. .\visualize.jl
```

A desktop with OpenGL support is required for the GLMakie window. Move the time slider or press Play. Close the window to exit. To check the mathematics without graphics packages:

```powershell
& $Julia --project=. .\test\runtests.jl
```

## What the visualization shows

- A reproducible random seed selects `(u0, v0)`; a small, nonzero square patch is present at time zero. The seed can be changed in `visualize.jl`.
- Time controls domain growth, not physical deformation: the surface stays on the same underlying torus. First the `u` range grows into a ribbon; at halfway, opposite `u` edges are identified. Then the `v` range grows; at the end its opposite edges are identified, producing a torus. Subdivision is represented by a fixed-resolution mesh sampled over the growing domain.
- The 2D panel shows the growing parameter rectangle. Orange and cyan mark the `u` and `v` edge pairs in both views. The parameter rectangle keeps both edges visible after identification; the surface shows a single seam instead. The status text reports when each identification has occurred.
- Red and blue line segments at the seed show scaled `x_u` and `x_v` tangent directions; they are not unit vectors.
- Surface color represents `|x_u × x_v|`, the local area scale per unit parameter area. The color range stays fixed throughout the animation.

## The builder

For major radius `R` and tube radius `r`, with `R > r > 0`, the reference map is

```text
x(u,v) = ((R + r cos(v)) cos(u), (R + r cos(v)) sin(u), r sin(v))
J(u,v) = [x_u  x_v]
|x_u × x_v| = r (R + r cos(v))
```

The map supplies the seed position; each other mesh position is built by integrating Jacobian columns along a path from the seed: first in `u`, then in `v`. Composite Simpson steps approximate these integrals; the test suite compares the result with the reference map and checks convergence. The Jacobian describes local changes, not the map's absolute position or its topology.

Closure is an extra topological rule: `(u,v) ~ (u+2pi,v)` and `(u,v) ~ (u,v+2pi)`. At each closure the mesh omits the duplicate endpoint and wraps triangle indices across the seam. The quadrature is approximate, so integrated values can differ slightly from exact periodic positions; the seam uses shared vertex indices, rather than pretending those approximations are identical. The initial parameter domain is a patch; the final quotient is a torus.

The spatial Jacobian has two columns. Time determines the parameter domain here; it is not included as a third Jacobian column or a fourth spatial dimension. A deforming map `x(u,v,t)` and its velocity `x_t` are possible later extensions.
