"""
Build and visualize an inversion model from perplex results. 

- plot pc1 vs comp w std of each bin 
- plot scatterplot of comp/seismic
"""

using DelimitedFiles
using HDF5
using StatGeochem
using ProgressMeter: @showprogress
using ArgParse
using Plots; gr();

include("../src/inversionModel.jl")
include("../src/config.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="test"
end
parsed_args = parse_args(ARGS, s)

models = makeModels(parsed_args["data_prefix"], modelType=InversionModel) 

@showprogress for layer in LAYER_NAMES
	for bin in 1:models.nbins
		model = models.models[layer][bin]
		p = scatter(model.seismic, model.comp, markeralpha=(min(1,1000/length(model.seismic)))); 
		outputPath = "data/"*parsed_args["data_prefix"]*"/output/model/"
		mkpath(outputPath) # make if does not exist
		savefig(p, "$outputPath/$(layer)_model_bin_$bin.png") # as a vector file it's too darn big
	end

	for bin in 1:models.nbins
		model = models.models[layer][bin]
		p = plot(model.c, model.m, yerror=model.e); 
		outputPath = "data/"*parsed_args["data_prefix"]*"/output/inversion_model/"
		mkpath(outputPath)
		savefig(p, "$outputPath/$(layer)_model_means_bin_$bin.pdf")
	end
end