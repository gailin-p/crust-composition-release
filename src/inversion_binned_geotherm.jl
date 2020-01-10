""" 
Inversion for binned geotherm 
(Eventually refactor for all inversions?)

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
        default="bin_geotherms"
end
parsed_args = parse_args(ARGS, s)

# Elements in ign files: index, elts, geotherm, layers. TODO add header to CSV created by resampleEarthChem
elements = ["index","SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2","padding", "tc1Crust","upper","middle","lower"] 

# Build models for every (geotherm bin)/layer combo 
filePrefix = "perplex_out_"
fileNames = filter(x->contains(x,"perplex_out_"), readdir("data/$(parsed_args["data_prefix"])"))
nBins = length(fileNames)
bins = crustDistribution.binBoundaries(nBins)

# Models for each layer. upper[i] is model for i-th geotherm bin. 
upper = Array{InversionModel, 1}(undef, nBins)
middle = Array{InversionModel, 1}(undef, nBins)
lower = Array{InversionModel, 1}(undef, nBins)

@showprogress 1 "Building models" for bin_num in 1:nBins
	# Build models for each layer for this bin 
	ignFile = "data/"*parsed_args["data_prefix"]*"/bsr_ignmajors_$(bin_num).csv"
	ign = readdlm(ignFile, ',')

	# Read perplex results. Shape (4, 3, n), property, layer, index. 
	perplexFile = "data/"*parsed_args["data_prefix"]*"/perplex_out_$(bin_num).h5"
	perplexresults = h5read(perplexFile, "results")

	if size(ign,1) != size(perplexresults,3)
		throw(AssertionError("Size of ign does not match size of perplex results from $(fileName)"))
	end 

	# Build models. use SiO2 (index 2 in ign). 
	upper[bin_num] = InversionModel(ign[:,1:2], Array(perplexresults[:,1,:]'))
	middle[bin_num] = InversionModel(ign[:,1:2], Array(perplexresults[:,2,:]'))
	lower[bin_num] = InversionModel(ign[:,1:2], Array(perplexresults[:,3,:]'))
end 

models = [upper, middle, lower]

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

	for (i, bin_bottom) in enumerate(bins[1:end-1])
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

	# save upper, middle, lower results 
	age_results[1,age_index] = nanmean(results[1])
	age_results[2,age_index] = nanmean(results[2])
	age_results[3,age_index] = nanmean(results[3])
end 


p = plot(size=(800,600));

layer_names = ["upper", "middle", "lower"]
for i in 1:3 # layers 
	plot!(p, ages, age_results[i,:], label=layer_names[i])
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




























