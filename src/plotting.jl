function h2plot(fom_name, method_vec, type="h2ham"; kwargs...)
    config = wload(datadir(fom_name, "config.jld2"))
    ylabel = type == "h2" ? L"Error $\|\Sigma_{\mathrm{pH}} - \widetilde{\Sigma}_{\mathrm{pH}}\|_{\mathcal{H}_2}$" : L"Error $\|\Sigma_{\mathcal{H}} - \widetilde{\Sigma}_{\mathcal{H}}\|_{\mathcal{H}_2}$"
    reduced_orders = config["reduced_orders"]
    p = plot(xticks=reduced_orders, xlimits=(reduced_orders[1],reduced_orders[end]), xlabel=L"Reduced order $r$", 
        ylabel=ylabel, yaxis=:log10, palette=:seaborn_bright; kwargs...)

    for method in method_vec
        h2errors_path = datadir(joinpath(fom_name, "h2errors"), savename(Dict("method" => method.name), "jld2"))
        try 
            h2errors = wload(h2errors_path)
            plot!(p, h2errors["r"], h2errors[type], label=method.label, marker=method.marker_shape, 
                color=method.color, markerstrokecolor=method.color, linestyle=method.linestyle)
        catch
            method_name = method.name
            @warn "Skipping plot for $method_name, since file $h2errors_path not exists."
            continue
        end
    end

    return p
end

function trajectoryplot(fom_name, method_vec, type="y", indexoutput=1; kwargs...)
    ylabels = Dict("y" => "Output", "ye" => "Output error", "yh" => "Hamiltonian", "yhe" => "Hamiltonian error")

    p = plot(xlabel="Time", ylabel=ylabels[type]; kwargs...)

    if type == "y" || type == "yh"
        fomsim = wload(datadir(joinpath(fom_name, "sims"), "fom.jld2"))
        plot!(p, fomsim["t"], fomsim[type][indexoutput,:], label="FOM", color=:black)
    elseif type == "ye" || type == "yhe"
        plot!(p, yaxis=:log10)
    end

    for method in method_vec
        sim = wload(datadir(joinpath(fom_name, "sims"), savename(Dict("method" => method.name), "jld2")))
        if type == "y" || type == "yh"
            plot!(p, sim["t"], sim[type][indexoutput, :], label=method.label, 
                color=method.color, linestyle=:dash)
        else
            plot!(p, sim["t"], sim[type][1, :], label=method.label, 
                color=method.color, linestyle=method.linestyle)
        end
    end

    return p
end
