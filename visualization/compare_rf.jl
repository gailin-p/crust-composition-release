"""
Run inversion on Crust1.0 seismic data and R&F seismic data 
"""

using DelimitedFiles
using HDF5
using StatGeochem
using ProgressMeter: @showprogress
using ArgParse
using Plots; gr();

include("../src/inversionModel.jl")
include("../src/invertData.jl")
include("../src/bin.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/bin_geotherm_si_weighted"
end
parsed_args = parse_args(ARGS, s)

function run(parsed_args)
	models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
	models = models.models # just interested in using individual models here 
	if length(models[UPPER]) > 1 
		throw(AssertionError("Can't use geotherm binning with R&F data"))
	end 

	# Get data to invert (resampled Crust1.0 data)
	(_, vp_upper, _, _),_ = getAllSeismic(6, n=2000)
	(_, vp_middle, _, _),_ = getAllSeismic(7, n=2000)
	(_, vp_lower, _, _),_ = getAllSeismic(8, n=2000)

	rf_upper = FRSeismic(6, n=2000)
	rf_middle = FRSeismic(7, n=2000)
	rf_lower = FRSeismic(8, n=2000)


	# Result data per geotherm bin 
	results_upper, errors_upper = estimateComposition(models[UPPER][1], vp_upper)
	results_middle, errors_middle  = estimateComposition(models[MIDDLE][1], vp_middle)
	results_lower, errors_lower = estimateComposition(models[LOWER][1], vp_lower)
	rf_results_upper, rf_errors_upper = estimateComposition(models[UPPER][1], rf_upper)
	rf_results_middle, rf_errors_middle  = estimateComposition(models[MIDDLE][1], rf_middle)
	rf_results_lower, rf_errors_lower = estimateComposition(models[LOWER][1], rf_lower)

	results = [results_upper, results_middle, results_lower]
	rf_results = [rf_results_upper, rf_results_middle, rf_results_lower]
	errors = [errors_upper, errors_middle, errors_lower]
	rf_errors = [rf_errors_upper, rf_errors_middle, rf_errors_lower]

	SI_index = findfirst(isequal("SiO2"), PERPLEX_ELEMENTS[2:end])
	println("Mean of upper results $(nanmean(results[1][:,SI_index]))")
	println("Mean of middle results $(nanmean(results[2][:,SI_index]))")
	println("Mean of lower results $(nanmean(results[3][:,SI_index]))")

	println("Mean of rf upper results $(nanmean(rf_results[1][:,SI_index]))")
	println("Mean of rf middle results $(nanmean(rf_results[2][:,SI_index]))")
	println("Mean of rf lower results $(nanmean(rf_results[3][:,SI_index]))")

	outputPath = "data/$(parsed_args["data_prefix"])/output"
	# Write to output files (one per layer)
	for (l, layer) in enumerate(LAYER_NAMES)
		output = hcat(results[l], errors[l])
		nelts = size(results[l],2)
		header = hcat(reshape(PERPLEX_ELEMENTS[2:nelts+1], (1,nelts)), # inverted data (results)
			reshape(PERPLEX_ELEMENTS[2:nelts+1], (1,nelts)).*" error") # errors on inverted data 
		output = vcat(header, output)
		writedlm(outputPath*"/results-rf_vs_crust1-crust1-$(layer).csv", output, ',')
	
		output = hcat(rf_results[l], rf_errors[l])
		nelts = size(results[l],2)
		header = hcat(reshape(PERPLEX_ELEMENTS[2:nelts+1], (1,nelts)), # inverted data (results)
			reshape(PERPLEX_ELEMENTS[2:nelts+1], (1,nelts)).*" error") # errors on inverted data 
		output = vcat(header, output)
		writedlm(outputPath*"/results-rf_vs_crust1-rf-$(layer).csv", output, ',')
	
	end 

end 

run(parsed_args)


























