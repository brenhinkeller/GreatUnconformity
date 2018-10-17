## -- Load external packages
    if VERSION>=v"0.7"
        using Statistics
        using DelimitedFiles
        using SpecialFunctions
    end

    using Plots; gr(); default(fmt=:svg);
    using StatGeochem

## --- Read Ronov volume data
    ronov = readdlm("RonovVolumes.csv",'|')
    ronov = elementify(ronov)

    x = vcat(ronov["t_age"]',ronov["b_age"]')[:]
    dt = ronov["b_age"]-ronov["t_age"]
    dt[1] = NaN # Exclude quaternary alluvium

    h = plot(x, repeat(ronov["Ronov_globe"]./dt,inner=2),label="Global")
    plot!(h, x, repeat(ronov["Ronov_NAm"]./dt,inner=2),label="N.Am.")
    plot!(h, x, repeat((ronov["Ronov_globe"] .- ronov["Ronov_NAm"])./dt,inner=2),label="Globe-N.Am.")
    plot!(h, xlims=(0,1600), xlabel="Age (Ma)", ylabel="Volume (km3/yr)", xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    savefig(h,"RonovVolume[Globe-NAm].pdf")
    display(h);

## --- Plot each continent by itself
    cont = ["Ronov_Afr", "Ronov_Eur", "Ronov_NAm", "Ronov_SAm", "Ronov_Aus"]
    h = plot(xlims=(0,1600),xlabel="Age (Ma)", ylims=(0,2),ylabel="Volume (km3/yr)")
    for i=1:5
        plot!(h, x, repeat(ronov[cont[i]]./dt,inner=2),label=continents[i], color=continentcolors[i])
    end
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    savefig(h,"RonovVolumebyContinent.pdf")
    display(h);

## --- Stacked histogram
    cont = ["Ronov_Afr", "Ronov_Eur", "Ronov_NAm", "Ronov_SAm", "Ronov_Aus"]
    vol_cumulative = Array{Float64}(undef, length(x),length(cont))
    total=fill(0.0, size(x))
    for i=5:-1:1
        total .+= repeat(ronov[cont[i]]./dt, inner=2)
        vol_cumulative[:,i] = total
    end

    h = plot()
    for i=1:5
        plot!(h,x,vol_cumulative[:,i],fill=0,label=continents[i],color=continentcolors[i])
    end
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    plot!(h, xlims=(0,1600), xlabel="Age (Ma)", ylims=(0,3.1), ylabel="Volume (km3/yr)")
    savefig(h,"RonovVolumebyContinentCumulative.pdf")
    display(h);

## --- Plot each continent by itself, normalized by area
    cont = ["Ronov_Afr", "Ronov_Eur", "Ronov_NAm", "Ronov_SAm", "Ronov_Aus"]
    contarea = [30.3E6, 53.4E6, 24.228E6, 18.28E6, 8.8015E6, 1.4E7]
    h = plot(xlims=(0,1600),xlabel="Age (Ma)", ylims=(0,4),ylabel="Volume flux (E-8 km3/yr/km2)")
    for i=1:5
        plot!(h, x, repeat(ronov[cont[i]]./dt/contarea[i]*1E8,inner=2), label=continents[i], color=continentcolors[i])
    end
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    savefig(h,"RonovVolumeFluxbyContinent.pdf")
    display(h);

## --- Compare Ronov and Macrostrat, for North America and globally

    # Load Macrostrat
    data = readdlm("../macrostrat/macrostrat.csv",',')
    macrostrat = elementify(data)

    # Compare Ronov and Macrostrat for North America
    h = plot(xlims=(0,1600),xlabel="Age (Ma)", ylims=(0,0.5),ylabel="Volume (km3/yr)")
    plot!(h, x, repeat(ronov["Ronov_NAm"]./dt,inner=2),label="Ronov (North America)", color=continentcolors[3])
    plot!(h, macrostrat["Time"], macrostrat["Volume"]/1E6,label="Macrostrat (North America)",color=continentcolors[3])
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)

    plot!(h,ylims=(0.006,5),yscale=:log)

    # Compare Ronov and Macrostrat scaled
    plot!(h, x, repeat(ronov["Ronov_globe"]./dt,inner=2),label="Ronov (global)", color=:black)
    plot!(h, macrostrat["Time"], macrostrat["Volume"]/1E6*sum(contarea)/contarea[3],label="Macrostrat scaled", color=:black)
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    display(h);

    savefig(h,"RonovVsMacrostrat.pdf")
    display(h);

## --- Plot glaciations
    h = plot(xlims=(0,1600),xflip=true,ylims=(0,1))
    plot!(h, [717,717],[0,1], label="")
    plot!(h, [635,635],[0,1], label="")
    plot!(h, [580,580],[0,1], label="")
    savefig(h, "glaciations.pdf")

## ---
