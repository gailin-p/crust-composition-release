""" 
Inversion for binned geotherm (works for not binned too; treats like nbins=1)

Main entry script for inversions. 
"""

using DelimitedFiles
using HDF5
using StatGeochem
using StatsBase
using ProgressMeter: @showprogress
using ArgParse
using Plots; gr();

include("inversionModel.jl")
include("invertData.jl")
include("../bin.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/latlong_weighted"
    "--model", "-m"
        help = "Type of model to use. Allowed: inversion (PCA), range (range of nearby samples)"
        arg_type = String
        range_tester = x -> (x in ["inversion","range"])
        default = "range"
    "--num_invert", "-n" 
    	help = "How many resampled Crust1.0 samples to invert?"
    	arg_type = Int
    	default = 50000
    "--age_model"
    	help = "How to apply ages to samples to invert. Allowed: tc1, earthchem"
        arg_type = String
        range_tester = x -> (x in ["tc1","earthchem"])
        default = "earthchem"
    "--bin_size", "-b"  
    	help = "% of each sizemic parameter in bin for range model (in decimal form)"
    	arg_type = Float64
    	default = .05
end
parsed_args = parse_args(ARGS, s)
writeOptions("data/"*parsed_args["data_prefix"]*"/inversion_options-$(parsed_args["model"])-$(parsed_args["age_model"]).csv", parsed_args)

function run(parsed_args)
	if parsed_args["model"] == "range"
		models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
		setError(models, parsed_args["bin_size"])
	elseif parsed_args["model"] == "inversion"
		models = makeModels(parsed_args["data_prefix"], modelType=InversionModel)
	end 

	# Get data to invert (resampled Crust1.0 data)
	# What age model to use? 
	if parsed_args["age_model"] == "tc1"
		age_model = Tc1Age()
	elseif parsed_args["age_model"] == "earthchem"
		age_model = EarthChemAge(10, 3)
	end

	upperDat, (upperCrustbase, upperAge, upperLat, upperLong) = getAllSeismic(6, n=parsed_args["num_invert"], ageModel=age_model, latlong=true) # returns rho, vp, vpvs, tc1, age
	middleDat, (middleCrustbase, middleAge, middleLat, middleLong) = getAllSeismic(7, n=parsed_args["num_invert"], ageModel=age_model, latlong=true)
	lowerDat, (lowerCrustbase, lowerAge, lowerLat, lowerLong) = getAllSeismic(8, n=parsed_args["num_invert"], ageModel=age_model, latlong=true)
	sampleAges = [upperAge, middleAge, lowerAge]
	sampleBases = [upperCrustbase, middleCrustbase, lowerCrustbase]
	sampleLats = [upperLat, middleLat, lowerLat]
	sampleLongs = [upperLong, middleLong, lowerLong]

	# Result data per geotherm bin 
	results_upper, errors_upper = estimateComposition(models, UPPER, upperDat...)
	results_middle, errors_middle  = estimateComposition(models, MIDDLE, middleDat...)
	results_lower, errors_lower = estimateComposition(models, LOWER, lowerDat...)
	results = [results_upper, results_middle, results_lower]
	errors = [errors_upper, errors_middle, errors_lower]

	SI_index = findfirst(isequal("SiO2"), PERPLEX_ELEMENTS[2:end])
	println("Mean of upper results $(nanmean(results[1][:,SI_index]))")
	println("Mean of middle results $(nanmean(results[2][:,SI_index]))")
	println("Mean of lower results $(nanmean(results[3][:,SI_index]))")

	# Oversampling ratio 
	n_original = length(crustDistribution.all_lats) # number of 1x1 grid cells w data at 
	n_resampled = parsed_args["num_invert"]

	outputPath = "data/"*parsed_args["data_prefix"]*"/output"
	mkpath(outputPath) # make if does not exist

	# Write to output files (one per layer)
	for (l, layer) in enumerate(LAYER_NAMES)
		#println("shapes $(size(sampleAges[l])), $(size(sampleBases[l])), $(size(results[l])), $(size(errors[l])),")
		output = hcat(sampleAges[l], sampleBases[l], sampleLats[l], sampleLongs[l], results[l], errors[l])
		nelts = size(results[l],2)
		header = hcat(["sample_age" "sample_depth" "sample_lat" "sample_long"], # sample data 
			reshape(PERPLEX_ELEMENTS[2:nelts+1], (1,nelts)), # inverted data (results)
			reshape(PERPLEX_ELEMENTS[2:nelts+1], (1,nelts)).*" error") # errors on inverted data 
		output = vcat(header, output)
		writedlm(outputPath*"/results-$layer-$(parsed_args["model"])-$(parsed_args["age_model"]).csv", output, ',')
	end 

end 

run(parsed_args)


























