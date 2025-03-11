using DrWatson
@quickactivate "ph-energy-matching"

using LinearAlgebra, Random, ControlSystems, ProgressMeter, Plots, LaTeXStrings
using PortHamiltonianBenchmarkSystems
using PortHamiltonianSystems, PortHamiltonianModelReduction, EnergyMatching, QuadraticOutputSystems

include(srcdir("experiments.jl"))
include(srcdir("methods.jl"))
include(srcdir("plotting.jl"))

import EnergyMatching: minreal

# set random seed
Random.seed!(1234)

# 1. Set up experiment
exp = msd(); # rcl(), msd(), poro(), msd_Xmin()
@info "Experiment: $(exp.name), n = $(ss(exp.fom).nx), m = $(ss(exp.fom).nu)"

# overwrites existing results in the data directory
overwrite = false
reduced_orders = 2:2:20

# 2. Apply minreal
Σ = minreal(exp.fom);
Σh = hdss(Σ);
h2 = norm(exp.fom - Σ)
h2ham = norm(hdss(exp.fom) - Σh)
@info "Minreal (n = $(ss(Σ).nx)): H2-error: $h2, Hamiltonian H2-error: $h2ham"
@tagsave(datadir(exp.name, "config.jld2"), Dict("exp" => exp.name, "reduced_orders" => reduced_orders, 
    "minreal" => Σ, "h2" => h2, "h2ham" => h2ham))


# 3. Run methods
reductors = [prbt_method, phirka_method]
matchers = [em_prbt_br_method, em_prbt_method, em_prbt_hypatia_method, em_phirka_method, em_phirka_hypatia_method]

# SDP solver which need extra installation and license, see https://github.com/jump-dev/MosekTools.jl
# You also need to uncomment the lines 5 and 27-28 in src/EnergyMatching/src/EnergyMatching.jl
# matchers = [em_prbt_mosek_method]

method_vec = vcat(reductors, matchers)

# Precompute positive real Gramians
Lx = prgrampd(Σ, :o)
norm(PortHamiltonianSystems.prare(Σ, :o, Lx'*Lx))/norm(ss(Σ).A)
Ly = prgrampd(Σ, :c)
norm(PortHamiltonianSystems.prare(Σ, :c, Ly'*Ly))/norm(ss(Σ).A)

for method in method_vec
    @info "Running " method.name

    @showprogress for (i,r) ∈ enumerate(reduced_orders)
        path = datadir(joinpath(exp.name, "roms"), 
            savename(Dict("method" => method.name, "r" => r), "jld2"))

        if !overwrite && isfile(path)
            @info "Skipping file already exists" path
            continue
        end

        if isa(method, Reductor)
            if method.name == "prbt"
                Σr = method.reduce(Σ, r; Lx=Lx, Ly=Ly)
            else
                Σr = method.reduce(Σ, r)
            end
        elseif isa(method, EnergyMatcher)
            rom_path = datadir(joinpath(exp.name, "roms"), savename(Dict("method" => method.rom_name, "r" => r), "jld2"))
            rom = wload(rom_path)["rom"]
            
            Σr = method.matchnrg(Σ, rom)
        else
            error("Method $method not recognized")
        end

        # save ROM
        @tagsave(path, Dict("rom" => Σr))
    end
end

# 4. Evaluate ROMs (H2-error and Hamiltonian H2-error)
for method in method_vec
    @info "Evaluating ROMs " method.name
    path = datadir(joinpath(exp.name, "h2errors"), savename(Dict("method" => method.name), "jld2"))
    
    if !overwrite && isfile(path)
        @info "Skipping file already exists" path
        continue
    end

    h2 = zeros(length(reduced_orders))
    h2ham = zeros(length(reduced_orders))
    @showprogress for (i,r) ∈ enumerate(reduced_orders) 
        rom_path = datadir(joinpath(exp.name, "roms"), savename(Dict("method" => method.name, "r" => r), "jld2"))
        rom = wload(rom_path)["rom"]
        
        h2[i] = norm(exp.fom - rom)
        h2ham[i] = h2norm(hdss(exp.fom) - hdss(rom))
    end
    @tagsave(path, Dict("method" => method.name, "r" => reduced_orders, "h2" => h2, "h2ham" => h2ham))
end

# 5. Analyzing results 

# H2-errors
method_vec_to_plot = [phirka_method, prbt_method, em_prbt_br_method, em_prbt_method, em_prbt_hypatia_method, em_phirka_method, em_phirka_hypatia_method]
p = h2plot(exp.name, method_vec_to_plot, "h2")
p = h2plot(exp.name, method_vec_to_plot, "h2ham")

# Trajectories

# input signal
u(x,t) = size(exp.fom.G, 2) == 1 ? [sin(t)] : [sin(t); cos(t)]
t = 0:0.1:100

res = lsim(ss(exp.fom), u, t; method=:RK4)
yh = sum(1//2 * res.x .* (exp.fom.Q * res.x); dims=1)
path = datadir(joinpath(exp.name, "sims"), "fom.jld2")
@tagsave(path, Dict("y" => res.y, "yh" => yh, "t" => res.t, "x" => res.x, "u" => res.u))

r = 16
for method in method_vec
    rom_path = datadir(joinpath(exp.name, "roms"), savename(Dict("method" => method.name, "r" => r), "jld2"))
    rom = wload(rom_path)["rom"]

    yr, _, xr, _ = lsim(ss(rom), u, t)
    yhr, _, _, _ = lsim(hdss(rom), u, t)
    # colwise norm
    ye = sqrt.(sum((res.y - yr).^2, dims=1))
    yhe = norm.(yh - yhr)

    path = datadir(joinpath(exp.name, "sims"), savename(Dict("method" => method.name), "jld2"))
    @tagsave(path, Dict("y" => yr, "yh" => yhr, "t" => t, "x" => xr, "u" => res.u, "ye" => ye, "yhe" => yhe))    
end

# plot trajectories
# ylims are set for msd, must be adjusted for other experiments
p = trajectoryplot(exp.name, method_vec, "y"; legend=:topright, xlims=(50,100), ylims=(-0.5,0.5))
p = trajectoryplot(exp.name, method_vec, "yh"; legend=:topright, xlims=(50,100), ylims=(0.2,0.7))
p = trajectoryplot(exp.name, method_vec, "ye"; legend=:topright, xlims=(50,100), ylims=(1e-8,1e-1))
p = trajectoryplot(exp.name, method_vec, "yhe"; legend=:topright, xlims=(50,100), ylims=(1e-2, 1e0))
