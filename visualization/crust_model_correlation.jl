"""
Invert for upper crust in specific 1x1 squares and plot against actual composition (from resampled ign)
for that 1x1 degree square. 

(To see if this model has high enough resolution/accuracy to differentiate between different crust squares)
"""

using DelimitedFiles 
using ArgParse
using Plots

include("../src/crustDistribution.jl")
include("../src/inversionModel.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        required= true
    "--model", "-m"
        help = "Type of model to use. Allowed: inversion, range"
        arg_type = String
        range_tester = x -> (x in ["inversion","range"])
        default = "inversion"
end
parsed_args = parse_args(ARGS, s)

# Build models 
if parsed_args["model"] == "range"
	models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
elseif parsed_args["model"] == "inversion"
	models = makeModels(parsed_args["data_prefix"], modelType=InversionModel)
end 

# Get crust distribution to invert 
rho, vp, vpvs, geotherms, lats, longs = crustDistribution.getAllSeismic(6, resample=true, latlong=true)

# Invert 
inverted_comp, inverted_errors = estimateComposition(models, UPPER, rho, vp, vpvs, geotherms)

# Get actual average comp for surface samples at each lat/long 
# Use resampled samples - TODO could consider re-doing resample, but seems unnececary  
ign = readdlm("data/$(parsed_args["data_prefix"])/bsr_ignmajors_1.csv", ',') # compositions same for all geotherm bins
lat_i = findfirst(isequal("Latitude"), PERPLEX_ELEMENTS)
long_i = findfirst(isequal("Longitude"), PERPLEX_ELEMENTS)
si_i = findfirst(isequal("SiO2"), PERPLEX_ELEMENTS)
surface_comp = fill(0.0, length(lats))
for i in 1:length(lats)
	lat = lats[i]
	long = longs[i]
	test = (abs.(ign[:,lat_i] .- lat) .< 1.0) .& (abs.(ign[:,long_i] .- long) .< 1.0)
	surface_comp[i] = mean(ign[test,si_i])
end 

p = scatter(surface_comp, inverted_comp, alpha=.3, legend=false, ylabel="Inverted SiO2 composition", xlabel="Surface SiO2 composition");
	#ylims=(30,90)); 


outputPath = "data/"*parsed_args["data_prefix"]*"/output"
savefig(p, "$outputPath/surface_vs_modeled-$(parsed_args["model"]).pdf")




