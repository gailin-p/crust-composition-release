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
using LoggingExtras

include("../src/config.jl")
include("../src/inversionModel.jl")
include("../src/invertData.jl")
include("../src/vpOnlyModel.jl")
include("../src/vpRhoModel.jl")
include("../src/linearModel.jl")
include("../src/rejectionModel.jl")
include("../src/rejectionModelSingle.jl")
include("../src/rejectionModelVpVpvs.jl")
include("../src/rejectionPriorModel.jl")
include("../src/bin.jl")
include("../src/rangeModel.jl")
include("../src/seismic.jl")
include("../src/parsePerplex.jl")

s = ArgParseSettings()
@add_arg_table! s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/base"
    "--name", "-o"
    	help = "Name of inversion run / folder to output results"
    	arg_type = String
    	default = "default"
    "--model", "-m"
        help = "Type of model to use. Allowed: inversion (PCA), range (range of nearby samples)"
        arg_type = String
        range_tester = x -> (x in ["inversion","range", "linear", "vprange", "vprhorange", "rejection", "vsrejection", "vprejection", "vpvsrejection", "priorrejection"])
        default = "rejection"
    "--num_invert", "-n"
    	help = "How many resampled Crust1.0 samples to invert?"
    	arg_type = Int
    	default = 252
    "--num_runs", "-r"
    	help = "How many times to run inversion?"
    	arg_type = Int
    	default = 1
    "--age_model"
    	help = "How to apply ages to samples to invert. Allowed: tc1, earthchem"
        arg_type = String
        range_tester = x -> (x in ["tc1","earthchem"])
        default = "earthchem"
    "--data_source"
    	help = "Source for seismic data"
        arg_type = String
        range_tester = x -> (x in ["Shen","Crust1.0", "Dabie", "DabieRG", "Test", "Spiral"])
        default = "Spiral"
	"--test_comp"
		help = "Composition of test earth"
		nargs = 3
		arg_type = Float64
		default = [66.6, 63.5, 53.4]
	"--test_source"
		help = "Data directory used as prior when building test earth. Leave blank to use inversion prior."
		arg_type = String
		default = ""
    "--data_source_uncertainty"
    	help = "Uncertainty as a fraction of std of this data set"
        arg_type = Float64
        range_tester = x -> (x >= 0)
        default = .1
    # "--mean"
    # 	help = "for range model, use the mean of all matching compsoitions? if not, choose best."
    # 	arg_type = Bool
    # 	default = false
    "--crack"
    	help = "Mean of upper crust total porosity."
    	arg_type = Float64
    	default = 0.007
    	range_tester = x -> ((x >= 0) && (x < .5))
    "--fraction_crack"
    	help = "% of the pore space that is cracks (very low aspect ratio).
    			Be careful with this, cracks have much larger effect on seismic props."
    	arg_type = Float64
    	default = 0.05
    	range_tester = x -> ((x >= 0) && (x <= 1))
    "--cracked_samples"
    	help = "How many samples do we apply cracking to?"
    	arg_type = Float64
    	default = 1.0
    	range_tester = x -> ((x >= 0) && (x <= 1))
	"--alteration_fraction"
		help = "Mean percent alteration"
		arg_type = Float64
		default = .01
		range_tester = x -> ((x >= 0) && (x <= .5))
end
parsed_args = parse_args(ARGS, s)
outputPath = "data/"*parsed_args["data_prefix"]*"/"*parsed_args["name"]*"/"
if ispath(outputPath)
	error("Output path $(outputPath) exists, delete or use a different --name option.")
end
if parsed_args["test_source"] == "" # default test data source is same as data source.
	parsed_args["test_source"] = parsed_args["data_prefix"]
end
mkpath(outputPath) # make output dir
writeOptions(outputPath*"/inversion_options.csv", parsed_args)

function run(parsed_args, outputPath)
	### SET UP CRACKING ###
	crackFile = outputPath*"crack_profile.csv"
	if !isfile(crackFile)
		println("Building new random cracking profiles...")
		ignFile = "data/"*parsed_args["data_prefix"]*"/bsr_ignmajors_1.csv"
		ign, header = readdlm(ignFile, ',', header=true)
		n = size(ign, 1)
		liquid_weights = [parsed_args["cracked_samples"]/2, parsed_args["cracked_samples"]/2, # dry, water
			0, 1- parsed_args["cracked_samples"]] # magma, no cracking
		crack_porosity = parsed_args["fraction_crack"] * parsed_args["crack"]
		pore_porosity = (1-parsed_args["fraction_crack"]) * parsed_args["crack"]

		if sum(liquid_weights) != 1
			error("Sum of crack weights does not equal 1. Is your cracked_samples > 1?")
		end

		profiles = Array{Crack,1}([
			random_cracking(crack_porosity, pore_porosity, parsed_args["alteration_fraction"],
			liquid_weights)
			for i in 1:n])
		write_profiles(profiles, crackFile)
	end

	### TODO: combine all these models, this is absurd -- most code is shared.
	if parsed_args["model"] == "range"
		models = makeModels(parsed_args["data_prefix"], modelType=RangeModel, crackFile=crackFile) # see inversionModel.jl. returns a ModelCollection
		# if parsed_args["mean"]
		# 	setMean(models, parsed_args["mean"])
		# end
	elseif parsed_args["model"] == "vprange"
		models = makeModels(parsed_args["data_prefix"], modelType=VpModel, crackFile=crackFile)
	elseif parsed_args["model"] == "vprhorange"
		models = makeModels(parsed_args["data_prefix"], modelType=VpRhoModel, crackFile=crackFile)
		# if parsed_args["mean"]
		# 	setMean(models, parsed_args["mean"])
		# end
	elseif parsed_args["model"] == "inversion"
		models = makeModels(parsed_args["data_prefix"], modelType=InversionModel, crackFile=crackFile)
	elseif parsed_args["model"] == "linear"
		models = makeModels(parsed_args["data_prefix"], modelType=LinearModel, crackFile=crackFile)
	elseif parsed_args["model"] == "rejection"
		models = makeModels(parsed_args["data_prefix"], modelType=RejectionModel, crackFile=crackFile)
	elseif parsed_args["model"] == "vprejection"
		models = makeModels(parsed_args["data_prefix"], modelType=VpRejectionModel, crackFile=crackFile)
	elseif parsed_args["model"] == "vsrejection"
		models = makeModels(parsed_args["data_prefix"], modelType=VsRejectionModel, crackFile=crackFile)
	elseif parsed_args["model"] == "vpvsrejection"
		models = makeModels(parsed_args["data_prefix"], modelType=RejectionModelVpVpvs, crackFile=crackFile)
	elseif parsed_args["model"] == "priorrejection"
		models = makeModels(parsed_args["data_prefix"], modelType=RejectionPriorModel, crackFile=crackFile)
	end

	# Get data to invert (resampled Crust1.0 data)
	# What age model to use?
	if parsed_args["age_model"] == "tc1"
		age_model = Tc1Age()
	elseif parsed_args["age_model"] == "earthchem"
		age_model = EarthChemAge(10, 3)
	end

	if parsed_args["data_source"] == "Test"
		source_models = makeModels(parsed_args["test_source"], # need cracked samples to build fake earths
				modelType=RejectionModel, crackFile=crackFile)
		upperDat, (upperCrustbase, upperAge, upperLat, upperLong) =
			getTestSeismic(parsed_args["num_invert"]*parsed_args["num_runs"], source_models, parsed_args["test_comp"][1], UPPER)
		middleDat, (middleCrustbase, middleAge, middleLat, middleLong) =
			getTestSeismic(parsed_args["num_invert"]*parsed_args["num_runs"], source_models, parsed_args["test_comp"][2], MIDDLE)
		lowerDat, (lowerCrustbase, lowerAge, lowerLat, lowerLong) =
			getTestSeismic(parsed_args["num_invert"]*parsed_args["num_runs"], models, parsed_args["test_comp"][3], LOWER)
	else
		upperDat, (upperCrustbase, upperAge, upperLat, upperLong) = getAllSeismic(
			6, n=parsed_args["num_invert"]*parsed_args["num_runs"], ageModel=age_model, latlong=true,
			dataSrc=parsed_args["data_source"], dataUncertainty=parsed_args["data_source_uncertainty"]) # returns rho, vp, vpvs, tc1, age
		middleDat, (middleCrustbase, middleAge, middleLat, middleLong) = getAllSeismic(
			7, n=parsed_args["num_invert"]*parsed_args["num_runs"], ageModel=age_model, latlong=true,
			dataSrc=parsed_args["data_source"], dataUncertainty=parsed_args["data_source_uncertainty"])
		lowerDat, (lowerCrustbase, lowerAge, lowerLat, lowerLong) = getAllSeismic(
			8, n=parsed_args["num_invert"]*parsed_args["num_runs"], ageModel=age_model, latlong=true,
			dataSrc=parsed_args["data_source"], dataUncertainty=parsed_args["data_source_uncertainty"])
	end
	sampleDat = [upperDat, middleDat, lowerDat]
	sampleAges = [upperAge, middleAge, lowerAge]
	sampleBases = [upperCrustbase, middleCrustbase, lowerCrustbase]
	sampleLats = [upperLat, middleLat, lowerLat]
	sampleLongs = [upperLong, middleLong, lowerLong]

	# Result data per geotherm bin
	results = [zeros(parsed_args["num_invert"]*parsed_args["num_runs"],length(PERPLEX_ELEMENTS)+1),
		zeros(parsed_args["num_invert"]*parsed_args["num_runs"],length(PERPLEX_ELEMENTS)+1),
		zeros(parsed_args["num_invert"]*parsed_args["num_runs"],length(PERPLEX_ELEMENTS)+1)]
	f_summary = outputPath*"runs_results.csv"

	## Run number of runs
	@showprogress "Running samples... " for i in 1:parsed_args["num_runs"]
		if parsed_args["model"] == "priorrejection"
			newPrior!(models)
		end

		istart = (i-1)*parsed_args["num_invert"] + 1
		iend = i*parsed_args["num_invert"]

		results_upper, _ = estimateComposition(models, UPPER,
			upperDat[1][istart:iend], upperDat[2][istart:iend], upperDat[3][istart:iend], upperDat[4][istart:iend])
		results_middle, _  = estimateComposition(models, MIDDLE,
			middleDat[1][istart:iend], middleDat[2][istart:iend], middleDat[3][istart:iend], middleDat[4][istart:iend])
		results_lower, _ = estimateComposition(models, LOWER,
			lowerDat[1][istart:iend], lowerDat[2][istart:iend], lowerDat[3][istart:iend], lowerDat[4][istart:iend])
		results[1][istart:iend,:] = results_upper
		results[2][istart:iend,:] = results_middle
		results[3][istart:iend,:] = results_lower

		SI_index = findfirst(isequal("SiO2"), PERPLEX_ELEMENTS)
		CO_index = findfirst(isequal("CO2"), PERPLEX_ELEMENTS)
		#println("Mean of upper results $(nanmean(results[1][:,SI_index:CO_index], dims=1))")
		#println("Mean of middle results $(nanmean(results[2][:,SI_index:CO_index], dims=1))")
		#println("Mean of lower results $(nanmean(results[3][:,SI_index:CO_index], dims=1))")

		io = open(f_summary, "a")
		writedlm(io, [nanmean(results[1][istart:iend,SI_index]),nanmean(results[2][istart:iend,SI_index]),nanmean(results[3][istart:iend,SI_index])]')
		close(io)
	end

	# Write to output files (one per layer)
	for (l, layer) in enumerate(LAYER_NAMES)
		#println("shapes $(size(sampleAges[l])), $(size(sampleBases[l])), $(size(results[l])), $(size(errors[l])),")
		output = hcat(sampleDat[l]...,
			sampleAges[l], sampleBases[l], sampleLats[l], sampleLongs[l], results[l])
		nelts = size(results[l],2) - 1 # last is geotherm bin info
		header = hcat(["sample_rho" "sample_vp" "sample_vpvs" "sample_geotherm"], # inverted data
			["sample_age" "sample_depth" "sample_lat" "sample_long"], # supplemental sample data
			reshape(PERPLEX_ELEMENTS, (1,nelts)), # inverted data (results)
			["bin"])

		if (parsed_args["crack"] > 0) & (parsed_args["model"] == "range") # only range model preserves 1:1 sample indexing
			(cracks, crack_header) = readdlm(crackFile, ',', header=true)
			out_cracks = Array{Any,2}(undef, (size(results[l],1), size(cracks,2)))
			for j in 1:size(out_cracks,1)
				idx = results[l][j,1]
				if isnan(idx)
					out_cracks[j,:] .= NaN
				else
					out_cracks[j,:] = cracks[Int(idx),:]
				end
			end
			output = hcat(output, out_cracks)
			header = hcat(header, reshape(crack_header, (1, length(crack_header))))
		end

		output = vcat(header, output)
		filename = outputPath*"results-$layer.csv"
		writedlm(filename, output, ',')
	end

	modelFile = outputPath*"model_info.csv"
	modelInfo = modelSummary(models)
	write(modelFile, modelInfo)

end

with_logger(FileLogger("$outputPath/crust_data.log")) do
	run(parsed_args, outputPath)
end
