"""
Simulations of the shallow water equations for the two-basin problem.

Matthew N. Crowe

Project with Michael Nguyen and Ted Johnson.

"""

using Plots, NetCDF

h =  ncread("shallow_water.nc", "h")[:, :, 1, :]

x = Array(range(-grid.Lx/2, grid.Lx/2, grid.Nx))
y = Array(range(-grid.Ly/2, grid.Ly/2, grid.Ny))

N = size(h)[3]

anim = @animate for i ∈ 1:N
    @info "frame = $(i)"
    heatmap(x, y, h[:, :, i] .- exp(-b*(σ₀ - σ₂));
        colormap = :balance,
        xlabel = "x",
        ylabel = "y",
        title = "h(x, y)",
        clims=(-0.15, 0.15),
        aspect_ratio = 1,
        xlim = (x[1], x[grid.Nx]),
        ylim = (y[1], y[grid.Ny]),
    )
end

mp4(anim, "h_mov.mp4"; fps = 15)
