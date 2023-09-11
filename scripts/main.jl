using DrWatson
@quickactivate "ph-energy-matching"

using LinearAlgebra, Random, ControlSystems, ColorBrewer, ProgressMeter
using Plots, LaTeXStrings, PGFPlotsX 
using PortHamiltonianBenchmarkSystems
using PortHamiltonianSystems, PortHamiltonianModelReduction, EnergyMatching, QuadraticOutputSystems

include(srcdir("models.jl"))
include(srcdir("methods.jl"))
include(srcdir("plotting.jl"))

# create directories if they don't exist
mkpath(datadir())
mkpath(plotsdir())

# set random seed
Random.seed!(1234)

# 1. Setup
# FOM
fom = msd() # msd(), poro(), msd_Xmin()

# overwrites existing results in the data directory
overwrite = true

reduced_orders = 2:2:20
@tagsave(datadir(fom.name, "config.jld2"), Dict("fom" => fom.name, "reduced_orders" => reduced_orders))

# methods
reductors = [phirka_method, prbt_method]
matchers = [em_prbt_br_method, em_prbt_bfgs_method, em_prbt_cosmo_method]

# SDP solver which need extra installation, see https://github.com/jump-dev/MosekTools.jl and https://github.com/jump-dev/SeDuMi.jl
# You also need to uncomment the lines 4 and 25-28 in src/EnergyMatching/src/EnergyMatching.jl
# matchers = [em_prbt_mosek_method, em_prbt_sedumi_method]

method_vec = vcat(reductors, matchers)

# precompute positive real Gramians 
# for the msd example the Riccati solution from MatrixEquations.jl works fine!
Lx = prgrampd(fom.Σ, :o)
Ly = prgrampd(fom.Σ, :c)

# for the poro example only the Riccati solution from the MATLAB solver icare works...
# for that you need to have MATLAB installed and the MATLAB.jl package,
# then uncomment line 4 in src/PortHamiltonianSystems/srcPortHamiltonianSystems.jl
# and uncomment the lines 72-77 in src/PortHamiltonianSystems/src/gramians.jl
# as well as the two lines below
# Lx = prgrampd(fom.Σ, :o; solver=:MATLAB) 
# Ly = prgrampd(fom.Σ, :c; solver=:MATLAB)

# 2. Run methods
for method in method_vec
    @info "Running " method.name

    @showprogress for (i,r) ∈ enumerate(reduced_orders)
        path = datadir(joinpath(fom.name, "roms"), 
            savename(Dict("method" => method.name, "r" => r), "jld2"))

        if !overwrite && isfile(path)
            @info "Skipping file already exists" path
            continue
        end

        if isa(method, Reductor)
            if method.name == "prbt"
                Σr = method.reduce(fom.Σ, r, Lx, Ly)
            else
                Σr = method.reduce(fom.Σ, r)
            end
        elseif isa(method, EnergyMatcher)
            rom_path = datadir(joinpath(fom.name, "roms"), savename(Dict("method" => "prbt", "r" => r), "jld2"))
            rom = wload(rom_path)["rom"]
           
            Σr = method.matchnrg(fom.Σ, rom)
        else
            error("Method $method not recognized")
        end

        # save ROM
        @tagsave(path, Dict("rom" => Σr))
    end
end

# 3. Evaluate ROMs
for method in method_vec
    @info "Evaluating ROMs " method.name
    path = datadir(joinpath(fom.name, "h2errors"), savename(Dict("method" => method.name), "jld2"))
    
    if !overwrite && isfile(path)
        @info "Skipping file already exists" path
        continue
    end

    h2 = zeros(length(reduced_orders))
    h2ham = zeros(length(reduced_orders))
    @showprogress for (i,r) ∈ enumerate(reduced_orders) 
        rom_path = datadir(joinpath(fom.name, "roms"), savename(Dict("method" => method.name, "r" => r), "jld2"))
        rom = wload(rom_path)["rom"]
        
        # Compute H2 norm and Hamiltonian H2 norm
        h2[i] = norm(fom.Σ - rom)
        h2ham[i] = norm(hdss(fom.Σ) - hdss(rom))
    end
    @tagsave(path, Dict("method" => method.name, "r" => reduced_orders, "h2" => h2, "h2ham" => h2ham))
end

# 4. Analyzing results
# H2 norms
p = h2_plot(fom.name, method_vec, "h2")
savefig(p, plotsdir(fom.name * "_h2.pdf"))
p = h2_plot(fom.name, method_vec, "h2ham")
savefig(p, plotsdir(fom.name * "_hamiltonian_h2.pdf"))

# trajectories
colwisenorm(x) = sqrt.(sum(x.^2, dims=1))

# input signal
u(x,t) = [sin(t); cos(t)]
# Time domain
t = 0:0.1:100

y, x, _, uout = lsim(ss(fom.Σ), u, t)
yh, _, _, _ = lsim(hdss(fom.Σ), u, t) 
path = datadir(joinpath(fom.name, "sims"), "fom.jld2")
@tagsave(path, Dict("y" => y, "yh" => yh, "t" => t, "x" => x, "u" => uout))

r = 20
for method in method_vec
    rom_path = datadir(joinpath(fom.name, "roms"), savename(Dict("method" => method.name, "r" => r), "jld2"))
    rom = wload(rom_path)["rom"]

    yr, _, xr, _ = lsim(ss(rom), u, t)
    yhr, _, _, _ = lsim(hdss(rom), u, t)
    ye = colwisenorm(y - yr)
    yhe = norm.(yh - yhr)

    path = datadir(joinpath(fom.name, "sims"), savename(Dict("method" => method.name), "jld2"))
    @tagsave(path, Dict("y" => yr, "yh" => yhr, "t" => t, "x" => xr, "uout" => u, "ye" => ye, "yhe" => yhe))    
end

# plot trajectories
# p = trajectory_plot(fom.name, method_vec, "y"; legend=:topright, xlimits=(50,100), ylimits=(-0.5,0.5))
# savefig(p, plotsdir(fom.name * "_output.pdf"))
p = trajectory_plot(fom.name, method_vec,  "yh"; legend=:topright, xlimits=(50,100))
savefig(p, plotsdir(fom.name * "_hamiltonian_output.pdf"))
p = trajectory_plot(fom.name, method_vec, "ye"; legend=:topright, xlimits=(50,100), ylimits=(1e-6,1e-1))
savefig(p, plotsdir(fom.name * "_output_error.pdf"))
p = trajectory_plot(fom.name, method_vec, "yhe"; legend=:topright, xlimits=(50,100), ylimits=(1e-3,1e0))
savefig(p, plotsdir(fom.name * "_hamiltonian_output_error.pdf"))