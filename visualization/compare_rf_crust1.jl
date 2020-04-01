"""
Compare vp distributions for Rudinck and Fountain and Crust1.0 
"""

using DelimitedFiles 
using Plots
gr();
using StatsPlots
include("../src/invertData.jl")

(rho, vp_upper, vpvs, geotherm),ages = getAllSeismic(6, n=1000)
(rho, vp_middle, vpvs, geotherm),ages = getAllSeismic(7, n=1000)
(rho, vp_lower, vpvs, geotherm),ages = getAllSeismic(8, n=1000)

rf_upper = FRSeismic(6)
rf_middle = FRSeismic(7)
rf_lower = FRSeismic(8)

# TODO grouped box plot? 
p = boxplot(["Upper"], [rf_upper], legend=false, color=:blue)
p = boxplot!(["Upper"], [vp_upper], color=:red)
p = boxplot!(["Middle"], [rf_middle], color=:blue)
p = boxplot!(["Middle"], [vp_middle], color=:blue)
p = boxplot!(["Middle"], [rf_middle], color=:red)
p = boxplot!(["Lower"], [rf_lower], color=:blue)
p = boxplot!(["Lower"], [vp_lower], color=:red)

outputPath = "data/output/RF_Crust-compare.pdf"
savefig(p, outputPath)