    using Plots; gr();
    using StatGeochem

    data = readcsv("GSCmaparea.csv")
    sedarea = elementify(data)

    h = plot(sedarea["age_Ma"],sedarea["sed_global_km2_Myr"],label="Global")
    plot!(h,sedarea["age_Ma"],sedarea["sed_Nam_km2_Myr"],label="N.Am.")
    plot!(h,sedarea["age_Ma"],sedarea["sed_non_Nam_km2_Myr"],label="Globe - N.Am.")
    plot!(h,xlims=(0,4000),ylims=(0,110000),xflip=true,legend=:topleft,fg_color_legend=:white,framestyle=:box,grid=false)
    plot!(h,xlabel="Age (Ma)", ylabel="Exposed bedrock area (km2/Myr)")
    savefig(h,"GSC Bedrock exposure, [meta]sedimentary.pdf")
    display(h)
    plot!(h,xlims=(0,2500))
    savefig(h,"GSC Bedrock exposure, [meta]sedimentary [0-2500].pdf")
    plot!(h,xlims=(0,1600))
    savefig(h,"GSC Bedrock exposure, [meta]sedimentary [0-1600].pdf")


## --- Read Ronov volume data

    ronov = readdlm("RonovVolumes.csv",'|')
    ronov = elementify(ronov)

    x = vcat(ronov["t_age"]',ronov["b_age"]')[:]
    dt = ronov["b_age"]-ronov["t_age"]
    dt[1] = NaN # Exclude quaternary alluvium

    h = plot(x, repeat(ronov["Ronov_globe"]./dt,inner=2),label="Global")
    plot!(h, x, repeat(ronov["Ronov_NAm"]./dt,inner=2),label="N.Am.")
    plot!(h, x, repeat((ronov["Ronov_globe"] .- ronov["Ronov_NAm"])./dt,inner=2),label="Globe-N.Am.")
    plot!(h, xlims=(0,1600), xlabel="Age (Ma)", ylabel="Volume (km3/Myr)", xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    savefig(h,"RonovVolume[Globe-NAm].pdf")
    display(h)

## --- Plot each continent by itself
    h = plot(x, repeat(ronov["Ronov_Eur"]./dt,inner=2),label="Eurasia")
    plot!(h, x, repeat(ronov["Ronov_NAm"]./dt,inner=2),label="N.Am.")
    plot!(h, x, repeat(ronov["Ronov_SAm"]./dt,inner=2),label="S.Am.")
    plot!(h, x, repeat(ronov["Ronov_Afr"]./dt,inner=2),label="Africa")
    plot!(h, x, repeat(ronov["Ronov_Aus"]./dt,inner=2),label="Australia")
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    plot!(h, xlims=(0,1600),xlabel="Age (Ma)", ylims=(0,2),ylabel="Volume (km3/Myr)")
    savefig(h,"RonovVolumebyContinent.pdf")
    display(h)

## --- Stacked histogram
    cont = ["Ronov_Aus","Ronov_SAm","Ronov_NAm","Ronov_Eur","Ronov_Afr",]
    vol_cumulative = Array{Float64}(length(x),length(cont))
    total=fill(0,size(x))
    for i=1:length(cont)
        total += repeat(ronov[cont[i]]./dt, inner=2)
        vol_cumulative[:,i] = total
    end

    h = plot()
    for i=length(cont):-1:1
        plot!(h,x,vol_cumulative[:,i],fill=0,label=cont[i])
    end
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    plot!(h, xlims=(0,1600), xlabel="Age (Ma)", ylims=(0,3.1), ylabel="Volume (km3/Myr)")
    savefig(h,"RonovVolumebyContinentCumulative.pdf")
    display(h)
## --- Plot each continent by itself, normalized by area

    h = plot(x, repeat(ronov["Ronov_Eur"]./dt/53.4E6*1E8,inner=2),label="Eurasia")
    plot!(h, x, repeat(ronov["Ronov_NAm"]./dt/24.228E6*1E8,inner=2),label="N.Am.")
    plot!(h, x, repeat(ronov["Ronov_SAm"]./dt/18.28E6*1E8,inner=2),label="S.Am.")
    plot!(h, x, repeat(ronov["Ronov_Afr"]./dt/30.3E6*1E8,inner=2),label="Africa")
    plot!(h, x, repeat(ronov["Ronov_Aus"]./dt/8.8015E6*1E8,inner=2),label="Australia")
    plot!(h, xflip=true, legend=:topleft, fg_color_legend=:white, framestyle=:box)
    plot!(h, xlims=(0,1600),xlabel="Age (Ma)", ylims=(0,4),ylabel="Volume flux (E-8 km3/Myr/km2)")
    savefig(h,"RonovVolumeFluxbyContinent.pdf")
    display(h)

## ---
