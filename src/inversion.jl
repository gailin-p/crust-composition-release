"""
Given a resampled EarthChem dataset and a perplex result file (in a data directory),
make an inversion model.

"""

using DelimitedFiles
using HDF5
using StatGeochem
using ProgressMeter: @showprogress
using ArgParse

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

# Build models. use SiO2. 
max_i = min(size(ign,1), size(perplexresults,3))
upper = InversionModel(ign[1:max_i,1:2], Array(perplexresults[:,1,1:max_i]'))
middle = InversionModel(ign[1:max_i,1:2], Array(perplexresults[:,2,1:max_i]'))
lower = InversionModel(ign[1:max_i,1:2], Array(perplexresults[:,3,1:max_i]'))


# Get seismic data, resample then invert 
# Upper 
rho, vp, vpvs = crustDistribution.getAllSeismic(6)

(means, errors) = estimateComposition(upper, rho, vp, vpvs)
println("Upper $(nanmean(means))")

# Middle 
rho, vp, vpvs = crustDistribution.getAllSeismic(7)
(means, errors) = estimateComposition(middle, rho, vp, vpvs)
println("Middle $(nanmean(means))")

# Lower
rho, vp, vpvs = crustDistribution.getAllSeismic(8)
(means, errors) = estimateComposition(middle, rho, vp, vpvs)
println("Lower $(nanmean(means))")


