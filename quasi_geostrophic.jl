"""
Simulations of the shallow water equations for the two-basin problem.

The basin shapes are set using the contours of the function

f(x, y) = |z² - 1|    where    z = x + iy,

This ensures that the problem can be transformed simply under the
conformal mapping transformation

σ + iτ = log(z² - 1).

The resulting (linearised) PV equation can now be solved for the
topographic Rossby wave modes in (σ, τ) space.

This code include the initial value problem with a Gaussian height
profile initial condition, and the forced problem with a Gaussian
forcing in the x velocity. These lines can be commented out as
required. Note that we solve the equations in conservative form so
our fields are h, uh and vh, rather than h, u and v.

Matthew N. Crowe

Project with Michael Nguyen and Ted Johnson.

"""

using GeophysicalFlows
using ImageFiltering
using NetCDF

# Set some physical and numerical parameters:

Nx, Ny = 1024, 1024                   # number of gridpoints
Lx, Ly = 4, 4                         # domain size
σ₀, σ₁, σ₂, b = 0.6, -0.5, -0.5, 1    # definitions for bottom bathymetry
H₀, L = 0.1, 0.1                      # height and width of initial perburbation
F₀ = 0                                # forcing magnitude
ν, nν = 0, 1                          # viscosity/hyperviscosity
μ = 0                                 # linear damping
δ = 0.01                              # smoothing lengthscale
T = 500                               # simulation stop time
saves = 1000                          # number of saves
dev = GPU()                           # run on CPU() or GPU()
Δt = 0.01                             # timestep
filename = "QG.nc"                    # output file

# Define number of iterations:

Nt = Int(T / Δt)
I = Int(Nt / saves)

# Define helper functions:

to_dev(f) = device_array(dev)(f)
to_CPU(f) = device_array(CPU())(f)

# Define the bottom topography:

grid = TwoDGrid(CPU(); nx=Nx, ny=Ny, Lx, Ly)    # build topography on CPU()
x, y = gridpoints(grid)

σ = @. log((x^2 - y^2 - 1)^2 + 4x^2*y^2) / 2

H = @. exp(-b*(σ₀ - σ₂)) * (σ >= σ₀) + (((σ >= σ₂) * exp(-b*(σ - σ₂)) + (σ < σ₂) * 1) * (x >= 0)
    + ((σ >= σ₁) * exp(-b*(σ - σ₂)) + (σ < σ₁) * exp(-b*(σ₁ - σ₂))) * (x < 0)) * (σ < σ₀)

N = Int(2* floor(δ / grid.dx) + 1)  # size of smoothing region
S = centered(ones(N, N) / N^2)      # smoothing filter

B = to_dev(imfilter(H, S, "circular"))

# Define PV forcing:

F = @. F₀ * exp(-(y^2 + (x-1.3)^2) / 0.2^2)    # physical space forcing
F_fft = to_dev(rfft(F))

function calcF!(Fh, sol, t, clock, vars, params, grid)
	@. Fh = F_fft
	return nothing
end

# Define problem:

prob = SingleLayerQG.Problem(dev;
                             nx=Nx,
                             ny=Ny,
                             Lx,
                             Ly,
			     dt=Δt,
                             stepper = "FilteredRK4",
                             aliased_fraction = 0,
                             calcF=calcF!,
                             eta=-B,
                             ν,
                             nν,
                             μ)

# Define and set initial condition:

q₀ = @. H₀ * exp(-((x-1.3)^2 + y^2)/L^2)

SingleLayerQG.set_q!(prob, to_dev(q₀))
SingleLayerQG.updatevars!(prob)

# Define diagnostics:

E = Diagnostic(SingleLayerQG.energy, prob; nsteps=Nt, freq=I)
Z = Diagnostic(SingleLayerQG.enstrophy, prob; nsteps=Nt, freq=I)
diags = [E, Z]

# Define output and save function:

if isfile(filename); rm(filename); end

t = LinRange(0, T, saves + 1)
nccreate(filename, "psi", "x", grid.x, "y", grid.y, "t", t)
nccreate(filename, "Q", "x", grid.x, "y", grid.y, "t", t)
nccreate(filename, "E", "t")
nccreate(filename, "Z", "t")

function Save_Data(filename, problem, diagnostics, i)

    ψ, q = reshape(to_CPU(problem.vars.ψ),(Nx, Ny, 1)), reshape(to_CPU(problem.vars.q),(Nx, Ny, 1))
    ncwrite(ψ, filename, "psi", start = [1, 1, i+1], count = [Nx, Ny, 1])
    ncwrite(q, filename, "Q", start = [1, 1, i+1], count = [Nx, Ny, 1])

    E, Z = diagnostics
    ncwrite([Lx*Ly*E.data[i+1]], filename, "E", start = [i+1], count = [1])
    ncwrite([Lx*Ly*Z.data[i+1]], filename, "Z", start = [i+1], count = [1])

    @info "Iteration: $(i * I), $(round(problem.clock.t / T * 100, sigdigits = 3))% complete, " *
    "t = $(round(problem.clock.t, sigdigits = 3))"

    return nothing

end

# Initial save:

Save_Data(filename, prob, diags, 0)

# Timestep problem:

for i in 1:saves

    stepforward!(prob, diags, I)
    SingleLayerQG.updatevars!(prob)

    Save_Data(filename, prob, diags, i)

    all(isfinite, prob.vars.q) ? nothing : (@error "NaN detected in fields")

end