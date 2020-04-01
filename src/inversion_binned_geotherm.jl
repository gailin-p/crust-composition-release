""" 
Inversion for binned geotherm (works for not binned too; treats like nbins=1)

"""

using DelimitedFiles
using HDF5
using StatGeochem
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
        default="remote/bin_geotherm_si_weighted"
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
    # "--bin_size", "-b" # TODO implement 
    # 	help = "% of each sizemic parameter in bin for range model (in decimal form)"
    # 	arg_type = Int
    # 	default = .05
end
parsed_args = parse_args(ARGS, s)
writeOptions("data/"*parsed_args["data_prefix"]*"/inversion_options-$(parsed_args["model"])-$(parsed_args["age_model"]).csv", parsed_args)

function run(parsed_args)
	if parsed_args["model"] == "range"
		models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
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
	n_original = length(crustDistribution.all_lats) # number of 1x1 grid cells w data at this age. 
	n_resampled = parsed_args["num_invert"]

	# By age! 

	ages = ageBins(age_model) 
	
	age_results = Array{Float64, 2}(undef, 3, length(ages)-1) # (layer, age)
	age_1std = Array{Float64, 2}(undef, 3, length(ages)-1) # (layer, age)

	for (age_index, age) in enumerate(ages[1:end-1]) # duplicate of above, but using age when requesting data to invert and saving mean data to age_results

		# Result data for this age bin 
		testUpper = (upperAge .>= age) .& (upperAge .< ages[age_index+1])
		testMiddle = (middleAge .>= age) .& (middleAge .< ages[age_index+1])
		testLower = (lowerAge .>= age) .& (lowerAge .< ages[age_index+1])

		# save upper, middle, lower results 
		age_results[1,age_index] = nanmean(results_upper[testUpper,SI_index]) 
		age_results[2,age_index] = nanmean(results_middle[testMiddle,SI_index])
		age_results[3,age_index] = nanmean(results_lower[testLower,SI_index])
		age_1std[1,age_index] = nanstd(results_upper[testUpper,SI_index]) * sqrt(n_original)/sqrt(n_resampled)
		age_1std[2,age_index] = nanstd(results_middle[testMiddle,SI_index]) * sqrt(n_original)/sqrt(n_resampled)
		age_1std[3,age_index] = nanstd(results_lower[testLower,SI_index]) * sqrt(n_original)/sqrt(n_resampled)
	end 


	p = plot(size=(800,600));

	for i in 1:3 # layers 
		plot!(p, ages[1:end-1], age_results[i,:], yerror=age_1std[i,:], label=LAYER_NAMES[i], markerstrokecolor=:auto)
	end 

	plot!(p, [200], [66.0], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:blue, label="Upper of R & F")
	plot!(p, [200], [60.6], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:orange, label="Middle of R & F")
	plot!(p, [200], [52.3], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:green, label="Lower of R & F")

	# Plot exposed 
	tmin = 0
	tmax = 4000
	nbins = 10
	samplePath = "data/$(parsed_args["data_prefix"])/bsr_ignmajors_1.csv"
	ign = readdlm(samplePath, ',')
	age_i = findfirst(isequal("Age"),PERPLEX_ELEMENTS)
	si_i = findfirst(isequal("SiO2"),PERPLEX_ELEMENTS)
	original = matread(IGN_FILE)
	age_centers, elt_means, elt_errors = bin(ign[:,age_i], ign[:,si_i],
	        tmin, tmax, length(ign[si_i])/length(original["SiO2"]), nbins)
	plot!(p, age_centers, elt_means, yerror=elt_errors, label="Exposed");

	plot!(p,xlabel="Age");
	plot!(p,ylabel="SiO2 composition");
	plot!(p,title="Composition estimations");
	outputPath = "data/"*parsed_args["data_prefix"]*"/output"
	mkpath(outputPath) # make if does not exist
	savefig(p, outputPath*"/composition-ages-$(parsed_args["model"])-$(parsed_args["age_model"]).pdf");


	# Plot composition by depth of crust base 
	minBase = minimum(minimum.(sampleBases))
	maxBase = maximum(maximum.(sampleBases))
	p = plot(xlabel="Base of Crust1.0 (km)", ylabel="Estimated silica (wt%)")
	for (l, layer) in enumerate(LAYER_NAMES)
		centers, aves, yerrors, _  = bin(sampleBases[l], results[l][:,SI_index], minBase, maxBase, 
			n_original/n_resampled, 10) # x,y,min,max,oversamplingratio,nbins
		plot!(p, centers, aves, yerror=yerrors, label=layer)
	end

	savefig(p, outputPath*"/composition-crustbase-$(parsed_args["model"])-$(parsed_args["age_model"]).pdf")


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


























