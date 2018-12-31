## --- Read Heuret traces

    # using Images, ImageView
    using FileIO: load
    using Images, Plots, Compat
    using StatGeochem

    img = load("TrenchSedimentThickness.tif")

    sedcolor = [RGB{N0f8}(0.208,0.325,0.643)
        RGB{N0f8}(0.29,0.608,0.835)
        RGB{N0f8}(0.063,0.722,0.925)
        RGB{N0f8}(0.388,0.78,0.78)
        RGB{N0f8}(0.631,0.831,0.651)
        RGB{N0f8}(0.776,0.878,0.639)
        RGB{N0f8}(0.976,0.929,0.51)
        RGB{N0f8}(1.0,0.792,0.392)
        RGB{N0f8}(0.98,0.667,0.29)
        RGB{N0f8}(0.949,0.4,0.404)
        RGB{N0f8}(0.937,0.224,0.369)
    ]

    # Trench sediment thickness from Heuret et al, 2012
    sedthick = [ 0.2
                 0.4
                 0.5
                 0.7
                 0.9
                 1.5
                 2.5
                 3.5
                 4.5
                 6.5
                 7.5]

    latimg = repmat(collect(linspace(-90,90,1125)),1,2250)
    lonimg = repmat(collect(linspace(-180,180,2250))',1125,1)

    t = img .== sedcolor[1]
    latlist = latimg[t]
    lonlist = lonimg[t]
    sedthicklist = ones(sum(t))*sedthick[1]
    for i=2:length(sedcolor)
        t = img .== sedcolor[i]
        latlist = vcat(latlist, latimg[t])
        lonlist = vcat(lonlist, lonimg[t])
        sedthicklist = vcat(sedthicklist, ones(sum(t))*sedthick[i])
    end

## --- Produce histogram of sed thickness vs trench length

    using StatsBase

    # Resample, to account uncertainty, weighting to account for change in length with latitude
    # Assume uncertainty is of the order of half the thickness
    trenchrepresentative = bsresample(sedthicklist,sedthicklist/8,100000,cos.(latlist))

    # Collect histogram
    binedges = 0:8
    trenchfilldist = fit(Histogram, trenchrepresentative[:,1], binedges).weights

    # Plot results
    trenchfilldist *= 65000/length(trenchrepresentative) # Normalize global trench length of about 65000 km
    h = plot(cntr(binedges),trenchfilldist,seriestype=:bar,bar_width=1,label="", framestyle=:box)
    plot!(h,xlims=(0,8),ylims=(0,4E4),xlabel="Trench sediment thickness (km)",ylabel="Trench length (km)", markerstrokecolor=:auto)
    savefig(h,"TrenchFillingExtent.pdf")
    display(h)

## --- Read eHf data with lat-lon constraints

    data = readdlm("Bataille2017S1.csv",'|')
    data = elementify(data)

    lat = data["latitude"]
    lon = data["longitude"]
    sed = Array{Float64}(size(lat))
    sed_sigma = Array{Float64}(size(lat))
    distance = Array{Float64}(size(lat))
    for i=1:length(lat)
        dist = arcdistance(lat[i],lon[i],latlist,lonlist)
        k = argmin(dist)
        sed[i] = sedthicklist[k]
        distance[i] = dist[k]
        # Assume uncertainty is of the order of half the thickness
        sed_sigma[i] = sedthicklist[k]/2
    end

    data["eHf"] = eHf(data["176Hf_177Hf_normalized"], data["176Lu_177Hf"], data["Age_UPb"])

## --- Correlate each zircon

    # Only look at samples younger than 100 Ma and within five arc degress of an arc
    t = (data["Age_UPb"] .< 100) .& (distance .< 5) .& (data["eHf"] .> -25)

    # Clean up data
    LuHf = data["176Lu_177Hf"][t]
    LuHfSigma = data["176Lu_177Hf±"][t]
    LuHfSigma[isnan.(LuHfSigma)] = nanmean(LuHfSigma)
    HfHf = data["176Hf_177Hf_normalized"][t]
    HfHfSigma = data["176Hf_177Hf±"][t]
    HfHfSigma[isnan.(HfHfSigma)] = nanmean(HfHfSigma)
    AgeUPb = data["Age_UPb"][t]
    Age_UPb_Sigma = 0.05*data["Age_UPb"][t]

    # Resample and plot
    # (c,m,el,eu) = bin_bsr_means(sed[t],data["eHf"][t],0,6,6,sed_sigma[t]/2,1000)
    (c,m,el,eu) = bin_bsr_eHf(sed[t],HfHf, LuHf, AgeUPb, 0,6,6, sed_sigma[t],HfHfSigma,LuHfSigma,Age_UPb_Sigma,10000)
    h = plot(c, m, yerror=(el,eu), seriestype=:scatter, label="", framestyle=:box, color=:darkred,markerstrokecolor=:auto)
    plot!(h, xlims=(0,6), xlabel="Present-day trench sediment thickness (km)", ylabel="Zircon epsilon Hf")
    savefig(h,"SedimentSubductionEHfCorrelation.pdf")
    display(h)

## ---
