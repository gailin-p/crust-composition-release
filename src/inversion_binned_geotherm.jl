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
end
parsed_args = parse_args(ARGS, s)

models = makeModels(parsed_args["data_prefix"]) # see inversionModel.jl. returns (upper, middle, lower) models.
nBins = length(models[1])
bins = crustDistribution.binBoundaries(nBins)

# Get data to invert (resampled Crust1.0 data)
upperDat = crustDistribution.getAllSeismic(6) # returns rho, vp, vpvs, tc1
middleDat = crustDistribution.getAllSeismic(7)
lowerDat = crustDistribution.getAllSeismic(8)
dat = [upperDat, middleDat, lowerDat]

# Result data per geotherm bin 
results_upper = Array{Float64, 1}(undef, length(upperDat[1]))
results_middle = Array{Float64, 1}(undef, length(middleDat[1]))
results_lower = Array{Float64, 1}(undef, length(lowerDat[1]))
results = [results_upper, results_middle, results_lower]

@showprogress "Running samples" for (i, bin_bottom) in enumerate(bins[1:end-1])
	# Look in ign for only the samples with geotherms in this bin 
	bin_top = bins[i+1]
	for layer in 1:3
		layerDat = dat[layer] 
		geotherm = layerDat[4]
		test = (geotherm .> bin_bottom) .& (geotherm .<= bin_top)
		(means, errors) = estimateComposition(models[layer][i], layerDat[1][test], layerDat[2][test], layerDat[3][test])
		results[layer][test] = means
	end 
end 

println("Mean of upper results $(nanmean(results[1]))")
println("Mean of middle results $(nanmean(results[2]))")
println("Mean of lower results $(nanmean(results[3]))")


# By age! 

ages = crustDistribution.all_ages[2:end,1] # first col is median age in bin, first row is NaN
age_results = Array{Float64, 2}(undef, 3, length(ages)) # (layer, age)
age_1std = Array{Float64, 2}(undef, 3, length(ages)) # (layer, age)

@showprogress 1 "running age samples:" for (age_index, age) in enumerate(ages) # duplicate of above, but using age when requesting data to invert and saving mean data to age_results

	# Get data to invert (resampled Crust1.0 data)
	upperDat = crustDistribution.getAllSeismic(6, age=age) # returns rho, vp, vpvs, tc1
	middleDat = crustDistribution.getAllSeismic(7, age=age)
	lowerDat = crustDistribution.getAllSeismic(8, age=age)
	dat = [upperDat, middleDat, lowerDat]

	# Result data per geotherm bin 
	results_upper = Array{Float64, 1}(undef, length(upperDat[1]))
	results_middle = Array{Float64, 1}(undef, length(middleDat[1]))
	results_lower = Array{Float64, 1}(undef, length(lowerDat[1]))
	results = [results_upper, results_middle, results_lower]
	error_upper = Array{Float64, 1}(undef, length(upperDat[1]))
	error_middle = Array{Float64, 1}(undef, length(middleDat[1]))
	error_lower = Array{Float64, 1}(undef, length(lowerDat[1]))
	inverted_error = [error_upper, error_middle, error_lower]

	for (i, bin_bottom) in enumerate(bins[1:end-1])
		# Look in ign for only the samples with geotherms in this bin 
		bin_top = bins[i+1]
		for layer in 1:3
			layerDat = dat[layer] 
			geotherm = layerDat[4]
			test = (geotherm .> bin_bottom) .& (geotherm .<= bin_top)
			(means, errors) = estimateComposition(models[layer][i], layerDat[1][test], layerDat[2][test], layerDat[3][test])
			results[layer][test] = means
			inverted_error[layer][test] = errors
		end 
	end 

	# save upper, middle, lower results 
	age_results[1,age_index] = nanmean(results[1])
	age_results[2,age_index] = nanmean(results[2])
	age_results[3,age_index] = nanmean(results[3])
	age_1std[1,age_index] = nanmean(inverted_error[1])
	age_1std[2,age_index] = nanmean(inverted_error[2])
	age_1std[3,age_index] = nanmean(inverted_error[3])
end 


p = plot(size=(800,600));

layer_names = ["upper", "middle", "lower"]
for i in 1:3 # layers 
	plot!(p, ages, age_results[i,:], yerror=age_1std[i,:], label=layer_names[i], markerstrokecolor=:auto)
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




























