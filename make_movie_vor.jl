"""
Simulations of the shallow water equations for the two-basin problem.

Matthew N. Crowe

Project with Michael Nguyen and Ted Johnson.

"""

using Plots, NetCDF

ω =  ncread("shallow_water.nc", "ω")[:, :, 1, :]

x = Array(range(-grid.Lx/2, grid.Lx/2, grid.Nx + 1))
y = Array(range(-grid.Ly/2, grid.Ly/2, grid.Ny + 1))

N = size(ω)[3]

anim = @animate for i ∈ 1:N
    @info "frame = $(i)"
    heatmap(x, y, ω[:, :, i] .- exp(-b*(σ₀ - σ₂));
        colormap = :balance,
        xlabel = "x",
        ylabel = "y",
        title = "ω(x, y)",
        clims=(-5, 5),
        aspect_ratio = 1,
        xlim = (x[1], x[grid.Nx]),
        ylim = (y[1], y[grid.Ny]),
    )
end

mp4(anim, "vor_mov.mp4"; fps = 15)
