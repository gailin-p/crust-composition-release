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
include("crustDistribution.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/bin_geotherm_si_weighted"
    "--model", "-m"
        help = "Type of model to use. Allowed: inversion, range"
        arg_type = String
        range_tester = x -> (x in ["inversion","range"])
        default = "inversion"
    "--num_invert", "-n"
    	help = "How many resampled Crust1.0 samples to invert?"
    	arg_type = Int 
    	default = 50000
end
parsed_args = parse_args(ARGS, s)

if parsed_args["model"] == "range"
	models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
elseif parsed_args["model"] == "inversion"
	models = makeModels(parsed_args["data_prefix"], modelType=InversionModel)
end 

# Get data to invert (resampled Crust1.0 data)
upperDat = crustDistribution.getAllSeismic(6, n=parsed_args["num_invert"]) # returns rho, vp, vpvs, tc1
middleDat = crustDistribution.getAllSeismic(7, n=parsed_args["num_invert"])
lowerDat = crustDistribution.getAllSeismic(8, n=parsed_args["num_invert"])

# Result data per geotherm bin 
results_upper, _ = estimateComposition(models, UPPER, upperDat[1], upperDat[2], upperDat[3], upperDat[4])
results_middle, _ = estimateComposition(models, MIDDLE, middleDat[1], middleDat[2], middleDat[3], middleDat[4])
results_lower, _ = estimateComposition(models, LOWER, lowerDat[1], lowerDat[2], lowerDat[3], lowerDat[4])
results = [results_upper, results_middle, results_lower]

println("Mean of upper results $(nanmean(results[1]))")
println("Mean of middle results $(nanmean(results[2]))")
println("Mean of lower results $(nanmean(results[3]))")


# By age! 

ages = crustDistribution.all_ages[2:end,1] # first col is median age in bin, first row is NaN
age_results = Array{Float64, 2}(undef, 3, length(ages)) # (layer, age)
age_1std = Array{Float64, 2}(undef, 3, length(ages)) # (layer, age)

@showprogress 1 "running age samples:" for (age_index, age) in enumerate(ages) # duplicate of above, but using age when requesting data to invert and saving mean data to age_results

	# Get data to invert (resampled Crust1.0 data)
	upperDat = crustDistribution.getAllSeismic(6, age=age, n=parsed_args["num_invert"]) # returns rho, vp, vpvs, tc1
	middleDat = crustDistribution.getAllSeismic(7, age=age, n=parsed_args["num_invert"])
	lowerDat = crustDistribution.getAllSeismic(8, age=age, n=parsed_args["num_invert"])

	# Result data for this age bin 
	results_upper, error_upper = estimateComposition(models, UPPER, upperDat[1], upperDat[2], upperDat[3], upperDat[4])
	results_middle, error_middle = estimateComposition(models, MIDDLE, middleDat[1], middleDat[2], middleDat[3], middleDat[4])
	results_lower, error_lower = estimateComposition(models, LOWER, lowerDat[1], lowerDat[2], lowerDat[3], lowerDat[4])

	# save upper, middle, lower results 
	age_results[1,age_index] = nanmean(results_upper)
	age_results[2,age_index] = nanmean(results_middle)
	age_results[3,age_index] = nanmean(results_lower)
	age_1std[1,age_index] = nanmean(error_upper) # TODO this is not right: std of a mean is not mean of the stds of samples
	age_1std[2,age_index] = nanmean(error_middle)
	age_1std[3,age_index] = nanmean(error_lower)
end 


p = plot(size=(800,600));

for i in 1:3 # layers 
	plot!(p, ages, age_results[i,:], yerror=age_1std[i,:], label=LAYER_NAMES[i], markerstrokecolor=:auto)
end 

plot!(p, [200], [66.0], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:blue, label="Upper of R & F")
plot!(p, [200], [60.6], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:orange, label="Middle of R & F")
plot!(p, [200], [52.3], seriestype=:scatter, marker=7,markershape=:star5, markercolor=:green, label="Lower of R & F")

plot!(p,xlabel="Age");
plot!(p,ylabel="SiO2 composition");
plot!(p,title="Composition estimations");
outputPath = "data/"*parsed_args["data_prefix"]*"/output"
mkpath(outputPath) # make if does not exist
savefig(p, outputPath*"/composition-ages-$(parsed_args["model"]).pdf");




























