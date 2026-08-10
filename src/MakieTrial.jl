using GLMakie
using LinearAlgebra

fig = Figure(figsize=(600, 600), fontsize=14)
ax = Axis(fig[1,1], )

arrows2d!(ax, [1], [1], [0, 1], [1, 1])