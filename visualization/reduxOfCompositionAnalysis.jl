## Using resampled dataset (w perplex data) from main analysis, redo figures

using Statistics, DelimitedFiles, MAT, StatsBase, Interpolations
using ProgressMeter: @showprogress
using Plots; gr();
using Logging
using MultivariateStats
using Random
using HDF5
using JLD
using ArgParse
using StatGeochem

include("../bin.jl")

# Import
println("importing")
mcign = Dict{String,Any}()
for name in ["Age", "Upper_PC1_SiO2", "Middle_PC1_SiO2", "Lower_PC1_SiO2", "SiO2",
    "Calc_Lower_Vp", "Calc_Lower_Rho", "Calc_Lower_VpVs",
    "Calc_Middle_Vp", "Calc_Middle_Rho", "Calc_Middle_VpVs",
    "Calc_Upper_Vp", "Calc_Upper_Rho", "Calc_Upper_VpVs",
    "Calc_Upper_PC1", "Calc_Lower_PC1", "Calc_Middle_PC1"]
    mcign[name] = h5read("../mcign_hpha.h5", name)
end

ign = matread("../igncn1.mat")

########## BELOW IS MOSTLY JUST COPIED OFF runPerplexCrustalStructure.jl
#######################
######################
########################

## --- Estimate compositions

iVar = "PC1"
elem = "SiO2"
nbins = 8
tmin = 0
tmax = 4000

# For each crust type
println("plotting")
layers=["Upper", "Middle", "Lower"]
#@showprogress "Running composition: " for crust in layers
p = plot(size=(800,600));
r = plot(size=(800,600), legend=:bottomright); # plot relationship between composition and PC1
for crust in layers

    # Plot relationship
    # Samples with non-null perplex results
    #element_test = map(!isnan, mcign["Calc_" * crust * "_" * iVar])
    element_test = map(!isnan, mcign["Calc_" * crust * "_Rho"] .+
            mcign["Calc_" * crust * "_Vp"] .+ mcign["Calc_" * crust * "_VpVs"])

    # Look at your independent (seismic) variable of choice
    x = mcign["Calc_" * crust * "_" * iVar]
    xmin = percentile(x, .5)
    xmax = percentile(x,99.5)

    # Bin your element of interest as a function of your seismic variable
    # center (seismic), mean (element of interest eg SiO2), error
    c, m, e, ydevs = bin(x, mcign[elem][element_test], xmin, xmax,
        length(mcign["SiO2"])/length(ign["SiO2"]), 40)

    # plot relationship for this age bin and layer
    plot!(r, c, m, yerror=e.*2, label=crust, markerstrokecolor=:auto);


    #% Plot the temporal trend for this crust type
    age_centers, elt_means, elt_errors, devs = bin(mcign["Age"], mcign[crust * "_" * iVar * "_" * elem],
            tmin, tmax, length(mcign["SiO2"])/length(ign["SiO2"]), nbins)
    plot!(p, age_centers, elt_means, yerror=elt_errors.*2, label=crust, markerstrokecolor=:auto);
end

# Plot exposed
age_centers, elt_means, elt_errors, devs = bin(mcign["Age"], mcign[elem],
        tmin, tmax, length(mcign["SiO2"])/length(ign["SiO2"]), nbins)
plot!(p, age_centers, elt_means, yerror=elt_errors.*2, label="Exposed", markerstrokecolor=:auto);

# Plot traditional upper-middle-lower values (66.0 - 60.6 - 52.3) from Table 9 of Rudnick and Fountain ’95

#hline!(p, [66.0], linestyle=:dot, linecolor=:blue, label="Upper of Rudnick and Fountain")
#hline!(p, [60.6], linestyle=:dot, linecolor=:red, label="Middle of Rudnick and Fountain")
#hline!(p, [52.3], linestyle=:dot, linecolor=:green, label="Lower of Rudnick and Fountain")

plot!(p, [200], [66.0], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:blue, label="Upper of Rudnick and Fountain")
plot!(p, [200], [60.6], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:orange, label="Middle of Rudnick and Fountain")
plot!(p, [200], [52.3], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:green, label="Lower of Rudnick and Fountain")

plot!(p,xlabel="Age");
plot!(p,ylabel="SiO2 composition");
plot!(p,title="Composition estimations");
savefig(p, "../output/composition-ages.pdf");

# Save relationship plot
plot!(r,xlabel="PC1 of perple_x seismic data");
plot!(r,ylabel="SiO2 composition");
plot!(r,title="PC1-composition relationships");
savefig(r, "../output/composition-pc1-relationship.pdf");


h = plot(size=(440,600),yticks=([1,2,3],["Lower","Middle","Upper"]), ylabel="Crust", xlabel="SiO2 [wt. %]", fg_color_legend=:white)
plot!(h,[57.2,67.5,58.7],[1,2,3], marker=5, label="This study", legend=:left)
plot!(h,[52.3,60.6,66.0],[1,2,3], marker=5, label="Traditional crustal structure\n(Rudnick & Fountain, 1995)")
plot!(h,[58.1],[3],seriestype=:scatter, marker=7,markershape=:star5,label="\nEarthChem exposed average")
plot!(h, xlims=(50, 70))
savefig(h,"../output/CrustalStructureComparison.pdf")
display(h)
