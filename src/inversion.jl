"""
Given a resampled EarthChem dataset and a perplex result file (in a data directory),
make an inversion model.

"""

using DelimitedFiles
using HDF5
using StatGeochem
using ProgressMeter: @showprogress
using ArgParse
using Plots; gr();

include("inversionModel.jl")
include("crustDistribution.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="test"
end
parsed_args = parse_args(ARGS, s)

# TODO move string constants etc into a utility file
perplexFile = "data/"*parsed_args["data_prefix"]*"/perplex_out.h5" 
ignFile = "data/"*parsed_args["data_prefix"]*"/bsr_ignmajors.csv"
prop_labels = ["rho,kg/m3","vp,km/s","vp/vs"] # properties in this order 

# Read perplex results. Shape (4, 3, n), property, layer, index. 
perplexresults = h5read(perplexFile, "results")

# Read ignmajors
ign = readdlm(ignFile, ',')
# index, elts, geotherm, layers. TODO add header to CSV created by resampleEarthChem
elements = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2","tc1Crust"]

# Build models. use SiO2 (index 2 in ign). 
max_i = min(size(ign,1), size(perplexresults,3))
upper = InversionModel(ign[1:max_i,1:2], Array(perplexresults[:,1,1:max_i]'))
middle = InversionModel(ign[1:max_i,1:2], Array(perplexresults[:,2,1:max_i]'))
lower = InversionModel(ign[1:max_i,1:2], Array(perplexresults[:,3,1:max_i]'))


# Get seismic data, resample then invert 
# Upper 
rho, vp, vpvs, unused = crustDistribution.getAllSeismic(6)

(means, errors) = estimateComposition(upper, rho, vp, vpvs)
println("Upper $(nanmean(means))")

# Middle 
rho, vp, vpvs, unused = crustDistribution.getAllSeismic(7)
(means, errors) = estimateComposition(middle, rho, vp, vpvs)
println("Middle $(nanmean(means))")

# Lower
rho, vp, vpvs, unused = crustDistribution.getAllSeismic(8)
(means, errors) = estimateComposition(lower, rho, vp, vpvs)
println("Lower $(nanmean(means))")



# Make plot of change over age 
layer_models = [upper, middle, lower]
layer_names = ["upper", "middle", "lower"]
crust_layers = [6,7,8]
ages = crustDistribution.all_ages[2:end,1] # first col is median age in bin, first row is NaN

p = plot(size=(800,600));
for (i, model) in enumerate(layer_models)
	comps = Array{Float64,1}(undef,length(ages))
	crust_layer = crust_layers[i]
	for (j, age) in enumerate(ages) 
		rho, vp, vpvs, unused = crustDistribution.getAllSeismic(crust_layer, age=age)
		(means, errors) = estimateComposition(layer_models[i], rho, vp, vpvs)
		comps[j] = nanmean(means)
	end
	plot!(p, ages, comps, label=layer_names[i])
end

plot!(p, [200], [66.0], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:blue, label="Upper of R & F")
plot!(p, [200], [60.6], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:orange, label="Middle of R & F")
plot!(p, [200], [52.3], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:green, label="Lower of R & F")

plot!(p,xlabel="Age");
plot!(p,ylabel="SiO2 composition");
plot!(p,title="Composition estimations");
outputPath = "data/"*parsed_args["data_prefix"]*"/output"
mkpath(outputPath) # make if does not exist
savefig(p, outputPath*"/composition-ages.pdf");









