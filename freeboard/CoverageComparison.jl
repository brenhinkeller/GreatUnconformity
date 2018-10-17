## --- Load external packages
    if VERSION>=v"0.7"
        using Statistics
        using DelimitedFiles
        using SpecialFunctions
    end

    using Plots; gr(); default(fmt=:svg);
    using StatGeochem

## --- Plot and compare different estimates of continental coverage

    # Area of the continents
    area_cont_km2 = 1.4894E8
    # Excluding antarctica
    area_cc_km2 = 1.4894E8 - 1.4E7

    # Digitized from scan of Egyed (1956) Nature with PlotDigitizer.
    # Paleozoic, Mesozoic, and "Caniozoic" rescaled to their modern boundaries (541,252,66,0 Myr)
    egyed_termier_myr = [529.076, 507.928, 486.289, 458.775, 440.531, 417.673, 394.023, 375.125, 361.855, 345.028, 330.946, 313.166, 306.193, 296.568, 282.748, 275.076, 258.413, 248.355, 227.428, 212.373, 198.061, 187.34, 175.453, 151.163, 139.336, 110.362, 85.9226, 73.1421, 64.2713, 52.789, 41.7133, 20.8795, 9.55939, 3.51455]
    egyed_termier_km2 = 1E6 .* [34.7499, 31.6007, 33.1326, 28.0706, 45.8783, 34.0062, 52.383, 55.2615, 28.7052, 42.7889, 32.8643, 32.09, 42.7662, 17.5917, 31.5325, 38.1309, 20.8666, 19.1131, 22.9933, 16.2493, 25.0294, 26.9874, 27.3472, 15.349, 29.2083, 36.0758, 34.9838, 30.7257, 15.583, 18.1485, 20.9982, 10.9498, 3.95325, 0.508505]

    egyed_strahov_myr = [484.608, 399.113, 329.952, 309.219, 279.509, 256.341, 225.505, 177.268, 151.933, 92.0874, 51.7846, 8.63077]
    egyed_strahov_km2 = 1E6 .* [65.7566, 72.6922, 58.6154, 53.227, 49.878, 19.2318, 26.9204, 39.7824, 24.1789, 54.748, 33.151, 1.23662]

    # From Ronov, 1992
    ronov_myr =  [3.9605, 14.1815, 28.465, 44.95, 61.0, 83.25, 122.75, 154.25, 168.8, 187.7, 219.15, 242.1, 249.685, 262.235, 285.6, 311.05, 341.05, 370.8, 388.0, 406.25, 423.3, 435.6, 451.1, 464.2, 477.7, 493.2, 506.0, 526.0]
    ronov_coverage = 1/100 .* [5, 10, 9, 18, 15, 29, 24, 24, 22, 18, 18, 18, 17, 18, 22, 25, 28, 33, 35, 28, 28, 38, 34, 39, 37, 33, 32, 30]

    # Macrostrat
    data = readdlm("../macrostrat/macrostrat.csv",',')
    macrostrat = elementify(data)

    # Plot results
    h = plot(xflip=true, xlims=(0,1600), ylims=(0, 0.5), fg_color_legend=:white, legend=:topleft, framestyle=:box)
    plot!(h, egyed_strahov_myr,egyed_strahov_km2/area_cont_km2,color=:lightblue,label="Egyed/Strahov (global)")
    plot!(h, egyed_strahov_myr,egyed_strahov_km2/area_cont_km2,color=:lightblue,label="", seriestype=:scatter, markersize=2, markerstrokecolor=:auto)

    plot!(h, egyed_termier_myr,egyed_termier_km2/area_cc_km2,color=:darkblue,label="Egyed/Termier (global)")
    plot!(h, egyed_termier_myr,egyed_termier_km2/area_cc_km2,color=:darkblue,label="", seriestype=:scatter, markersize=2, markerstrokecolor=:auto)

    plot!(h, ronov_myr, ronov_coverage, color=:black, label = "Ronov (global)")
    plot!(h, ronov_myr, ronov_coverage, color=:black, label = "", seriestype=:scatter,  markersize=2, markerstrokecolor=:auto)

    plot!(h, macrostrat["Time"], macrostrat["Coverage"], color=:darkred, label = "Macrostrat (N. Am.)")
    plot!(h, xlabel="Age (Ma)", ylabel="Fraction covered")
    savefig(h, "CoverageFraction.pdf")
    display(h);
