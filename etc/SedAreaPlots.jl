## -- Load external packages
    if VERSION>=v"0.7"
        using Statistics
        using DelimitedFiles
        using SpecialFunctions
    end

    using Plots; gr(); default(fmt=:svg);
    using StatGeochem

## --- Plot GSC map area
    data = readdlm("GSCmaparea.csv",',')
    sedarea = elementify(data)

    h = plot(sedarea["age_Ma"],sedarea["sed_global_km2_Myr"],label="Global",color=:darkred)
    plot!(h,sedarea["age_Ma"],sedarea["sed_Nam_km2_Myr"],label="N.Am.",color=continentcolors[3])
    plot!(h,sedarea["age_Ma"],sedarea["sed_non_Nam_km2_Myr"],label="Globe - N.Am.",color=:black)
    plot!(h,xlims=(0,4000),ylims=(0,110000),xflip=true,legend=:topleft,fg_color_legend=:white,framestyle=:box,grid=false)
    plot!(h,xlabel="Age (Ma)", ylabel="Exposed bedrock area (km2/Myr)")
    savefig(h,"GSC Bedrock exposure, [meta]sedimentary.pdf")
    display(h)
    plot!(h,xlims=(0,2500))
    savefig(h,"GSC Bedrock exposure, [meta]sedimentary [0-2500].pdf")
    plot!(h,xlims=(0,1600))
    savefig(h,"GSC Bedrock exposure, [meta]sedimentary [0-1600].pdf")

## ---
