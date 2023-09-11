function Plots.savefig(p, path; kwargs...)
    for (key, value) in kwargs
        if value === nothing
            continue
        end
    
        key = replace(String(key), "_" => " ")
        p.o.the_plot.elements[1].elements[1].options[key] = value
    end
    return pgfsave(path, p.o.the_plot)
end

function h2_plot(fom_name, methods, type="h2ham"; kwargs...)
    config = wload(datadir(fom_name, "config.jld2"))
    # ylabel = type == "h2" ? L"$\mathcal{H}_2$ error" : L"Energy $\mathcal{H}_{2}$ error"
    ylabel = type == "h2" ? L"Error $\|\Sigma_{\mathrm{pH}} - \widetilde{\Sigma}_{\mathrm{pH}}\|_{\mathcal{H}_2}$" : L"Error $\|\Sigma_{\mathcal{H}} - \widetilde{\Sigma}_{\mathcal{H}}\|_{\mathcal{H}_2}$"
    reduced_orders = config["reduced_orders"]
    p = plot(xticks=reduced_orders, xlimits=(reduced_orders[1],reduced_orders[end]), xlabel=L"Reduced order $r$", 
        ylabel=ylabel, yaxis=:log10; kwargs...)

    for method in methods
        h2errors_path = datadir(joinpath(fom_name, "h2errors"), savename(Dict("method" => method.name), "jld2"))
        try 
            h2errors = wload(h2errors_path)
            plot!(p, h2errors["r"], h2errors[type], label=latexstring(method.label), marker=method.marker_shape, 
                color=method.color, markerstrokecolor=method.color, linestyle=method.linestyle)
        catch
            method_name = method.name
            @warn "Skipping plot for $method_name, since file $h2errors_path not exists."
            continue
        end
    end

    return p
end

function trajectory_plot(fom_name, methods, type="y"; kwargs...)
    ylabel = nothing
    yaxis = nothing

    fom_data = wload(datadir(joinpath(fom_name, "sims"), "fom.jld2"))

    if type == "y"
        ylabel = L"Output $y$"
    elseif type == "yh"
        ylabel = L"Hamiltonian $\mathcal{H}$"
    elseif type == "ye"
        ylabel = L"Error $\|y - \widetilde{y}\|_2$"
        yaxis=:log10
    elseif type == "yhe"
        ylabel = L"Error $|y_\mathcal{H} - \widetilde{y}_{\mathcal{H}}|$"
        yaxis=:log10
    end

    p = plot(xlabel=L"Time $t$", ylabel=ylabel, yaxis=yaxis; kwargs...)

    if type == "y" || type == "yh"
        plot!(p, fom_data["t"], fom_data[type][1,:], label="FOM", 
            color=:black, linestyle=:solid, linewidth=1, thickness_scaling = 1)
    end
   
    for (i, method) in enumerate(methods)
        rom_data_path = datadir(joinpath(fom_name, "sims"), savename(Dict("method" => method.name), "jld2"))
        try 
            rom_data = wload(rom_data_path)
            # linsestyle = type == "y" || type == "yh" ? :dash : :solid
                
            plot!(p, rom_data["t"], rom_data[type]', label=latexstring(method.label),
                color=method.color, linestyle=method.linestyle, linewidth=1, thickness_scaling = 1)
        catch
            method_name = method.name
            @warn "Skipping plot for $method_name, since file $rom_data_path not exists."
            continue
        end
    end

    return p
end


