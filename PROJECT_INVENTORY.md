# JacobianTorus project inventory

A snapshot of the project layout and Julia environment. Byte counts are exact for the listed workspace files at inventory time. The project-file total excludes this report to avoid a self-referential report-size count.

## 1. Source inventory

The ten project files total **85,140 bytes**. Of that, implementation and tests are **11,907 bytes**; the rest is project configuration/lockfile and README documentation. This makes clear that 85 KB is not code alone.

```text
JacobianTorus/
|-- JacobianTorus.slnx                     210 bytes
|-- JacobianTorus.vcxproj                7,114 bytes
|-- JacobianTorus.vcxproj.filters        1,080 bytes
|-- JacobianTorus.vcxproj.user             168 bytes
|-- Manifest.toml                       61,195 bytes
|-- Project.toml                           185 bytes
|-- README.md                            3,281 bytes
|-- visualize.jl                         5,718 bytes
|-- src/
|   `-- TorusGeometry.jl                 3,982 bytes
`-- test/
	`-- runtests.jl                      2,207 bytes
```

- **Mathematical implementation:** `src/TorusGeometry.jl` (3,982 bytes).
- **Interactive visualization:** `visualize.jl` (5,718 bytes).
- **Regression tests:** `test/runtests.jl` (2,207 bytes; 44 checks pass).
- **Project documentation:** `README.md` (3,281 bytes).
- **Julia environment files:** `Project.toml` and `Manifest.toml`.
- **Visual Studio shell:** solution/project metadata exists; README describes the C++ project as empty, while the working prototype is Julia.

## 2. Generated artifacts

The local `.vs/` directory is generated IDE/cache state, not project source. At this snapshot it contained **41 files / 1,316,788 bytes** (about 1.26 MiB). No `bin/` or `obj/` build-output directories were present in the inventory. Including `.vs/`, the ten project files and cache total **1,401,928 bytes** (about 1.34 MiB); this excludes this report.

## 3. Julia environment

- **Julia used:** 1.12.0; project compatibility declares Julia 1.10 or later.
- **Manifest format:** 2.0, with **259 dependency entries**.
- **Direct package dependencies:** GLMakie 0.13.17 and GeometryBasics 0.5.13.

## Direct dependency map

```text
Project.toml
|-- GLMakie 0.13.17
|   |-- Makie 0.24.17                 # plotting, axes, heatmap, mesh, UI
|   |-- GLFW 4.0.0                    # windowing
|   |-- GeometryBasics 0.5.13         # points, meshes, triangle faces
|   |-- Observables 0.5.5             # reactive visualization state
|   |-- ShaderAbstractions 0.5.1      # rendering support
|   |-- Colors 0.13.4 / ColorTypes 0.12.5
|   `-- StaticArrays 1.9.25
`-- GeometryBasics 0.5.13
	|-- EarCut_jll 2.2.4+0
	|-- LinearAlgebra (Julia stdlib)
	|-- PrecompileTools 1.4.0
	|-- Random (Julia stdlib)
	`-- StaticArrays 1.9.25
```

The two direct dependencies share GeometryBasics; the manifest records their transitive packages and exact resolved versions. The tree above is intentionally abbreviated to the direct dependencies and notable immediate rendering/geometry dependencies rather than all 259 manifest entries.

## 4. Mathematical components

- **Torus parameterization:** `torus_point` maps the two angular parameters to a ring torus with major radius `R` and tube radius `r`.
- **Jacobian columns:** `tangents` computes `x_u` and `x_v`, the local directions spanning the surface tangent plane.
- **Local area scale:** `area_scale` computes `|x_u × x_v| = r(R + r cos(v))`.
- **Patch integration:** composite Simpson steps integrate the Jacobian columns from a seed point, first along `u`, then along `v`.
- **Visualization:** the parameter-domain heatmap and 3D surface use the same fixed color scale for the local area scale, linking parameter values to geometric area.

## 5. Topological components

- **Growth states:** the parameter patch expands in `u`, then `v`; each coordinate direction becomes closed at its corresponding stage.
- **Seam identification:** a closed direction omits its duplicate endpoint and wraps triangle indices across the seam.
- **Boundary behavior:** open directions expose two mesh boundaries; closed directions join the paired edges.
- **Result:** identifying both periodic parameter directions produces the torus topology. The seam rule is distinct from the numerical integration, whose quadrature can leave small positional differences at periodic endpoints.

## Useful commands

Run from the project root in PowerShell:

```powershell
$Julia = 'C:\Users\tmcla\AppData\Local\Programs\Julia-1.12.0\bin\julia.exe'
& $Julia --project=. .\test\runtests.jl
& $Julia --project=. .\visualize.jl
```
