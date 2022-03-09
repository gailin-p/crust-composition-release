"""
Define abstract model type and functions to build collections of models across geotherm bins.
"""

using MultivariateStats
using StatGeochem
using StatsBase
using Interpolations
using DelimitedFiles
using HDF5
using MAT
using Random
using Plots; gr();

include("bin.jl") # TODO use Brenhin's version (bin_bsr with nresamples=0)
include("crustDistribution.jl")
include("utilities.jl")
include("cracks.jl")

prop_labels = ["rho,kg/m3","vp,km/s","vp/vs"] # properties in this order

"""
Allow plug-in of different model types.
Every model must implement a constructor that takes ign and perplex data,
and an estimateComposition method that takes test seismic data, so that
any model can be interchangably used in a ModelCollection.
"""
abstract type AbstractModel end

"""
    InversionModel
Store information for inverting observed data according to perplex data for samples
Includes: PCA, mean/std for normalizing inputs, relationship between PCA and comp
Different models for each inversion (upper != middle != lower) and geotherm bin.
"""
struct  InversionModel <: AbstractModel
	# PCA of seismic vars
	rhoNorm::Norm
	vpNorm::Norm
	vpvsNorm::Norm
	pca::PCA
	# Binned relationship between PCA and composition
	interp::AbstractInterpolation # Interpolation from PCA to element of interest
	error_interp::AbstractInterpolation # Interpolation from PCA to standard deviation of elt of interest
	# The bins themselves
	c::AbstractRange
	m::Array{Float64,1}
	e::Array{Float64,1}
	# Original data
	seismic::Array{Float64,1} # Seismic data of perplex samples post-pca. length n
	comp::Array{Float64,1} # composition of element of interest. length n.
end

"""
	BadModel
A place-holder model for when we didn't have enough data
"""
struct BadModel <: AbstractModel end

"""
All the models for a dataset/layer combo. (so one of these for each upper/middle/lower)
If binned, models will contain nbins models.
"""
struct ModelCollection
	models::Dict{String,Array{AbstractModel,1}}
	bins::AbstractRange
	nbins::Integer
end

"""
Output any relevant model info. Assume any model within models will return same info.
"""
function modelSummary(models::ModelCollection)
	return modelSummary(models.models[UPPER][1])
end

"""
Output relevant model info. Default: empty string
"""
function modelSummary(model::AbstractModel)
	return "This model type did not implement model summary"
end

"""
	Explicitly set the error ranges for an InversionModel (not implemented )
"""
function setError(model::InversionModel, errors::Tuple{Float64,Float64,Float64})
	throw("bin size or margin setting not implemented for inversion model.")
end


"""
	Explicitly set the error ranges for an InversionModel (not implemented )
"""
function setMean(model::InversionModel, usemean::Bool)
	throw("mean setting does not exist for inversion model ")
end

"""
	Set whether to use mean in every model in a model collection of range models
"""
function setMean(models::ModelCollection, usemean::Bool)
	for layer in keys(models.models)
		for model in models.models[layer]
			setMean(model, usemean)
		end
	end
end

function resultSize(model::BadModel)
	return 1
end


"""
Composition estimate for a BadModel returns NaNs.
"""
function estimateComposition(model::BadModel,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})
	return (fill(NaN, length(rho)), fill(NaN, length(rho)))
end

"""
	estimateComposition

Estimate composition using a collection of geotherm-binned models.
"""
function estimateComposition(models::ModelCollection, layer::String,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1}, geotherms::Array{Float64, 1})

	if (length(rho) != length(vp)) | (length(vp) != length(vpvs))
		throw(AssertionError("Lengths of seismic datasets don't match"))
	end

	if !(layer in LAYER_NAMES)
		throw(AssertionError("Invalid layer, use one of $LAYER_NAMES"))
	end

	if (sum(isnan.(geotherms))>0)
		throw(AssertionError("NaN geotherm in inversion"))
	end

	# Input and result data for all geotherm bins
	nresults = resultSize(models.models["upper"][1])
	results = fill(0.0, (length(rho), nresults+1))
	errors = fill(0.0, (length(rho), nresults))
	done = fill(0, length(rho))

	for (i, bin_bottom) in enumerate(models.bins[1:end-1])
		# Look in ign for only the samples with geotherms in this bin
		bin_top = models.bins[i+1]

		if length(models.bins) == 2 # only one bin, run all geotherms
			test = fill(true, length(geotherms))
		elseif i == 1 # first bin; include all low geotherms
			test = (geotherms .<= bin_top)
		elseif i == (length(models.bins) - 1) # top bin, include all high geotherms
			test = (geotherms .> bin_bottom)
		else
			test = (geotherms .> bin_bottom) .& (geotherms .<= bin_top)
		end

		if sum(test) == 0
			continue # This bin is empty
		end
		(means, es) = estimateComposition(models.models[layer][i], rho[test], vp[test], vpvs[test])
		results[test,1:end-1] = means
		results[test,end] .= i # record bin
		errors[test,:] = es
		done[test] .= 1
	end

	if sum(done) != length(done)
		throw(AssertionError("Did not process all samples. Processed $(sum(done)) out of $(length(done))"))
	end

	return (results, errors)
end

"""
	makeModels(data_location)

Load data for an inversion run (function expects ign csv from resampleEarthChem and subsequent perplex results)

Returns touple of lists of models, (upper, middle, lower), each with as many models as geotherm bins.
"""
function makeModels(data_location::String; modelType::DataType=RangeModel, crackProfile::CrackProfile)
	# Be flexible about where we're running from
	if isdir("data")
		folder_prefix = "data"
	else
		folder_prefix = "../data"
	end

	# Elements in ign files: index, elts, geotherm, layers. TODO add header to CSV created by resampleEarthChem
	elements = PERPLEX_ELEMENTS

	# Build models for every (geotherm bin)/layer combo
	fileNames = filter(x->contains(x,"perplex_out_"), readdir("$folder_prefix/$data_location"))
	nBins = length(fileNames)
	if nBins == 0
		throw(AssertionError("No perplex results found for data directory $data_location"))
	end
	bins = crustDistribution.binBoundaries(nBins) # bin by depth to 550
	bin_centers = [(bins[b]+bins[b+1])/2 for b in 1:length(bins)-1]
	bin_geotherms = 550 ./ bin_centers
	bin_crack_adjustments = izvestiya_conversion(bin_geotherms)

	# Models for each layer. upper[i] is model for i-th geotherm bin.
	upper = Array{AbstractModel, 1}(undef, nBins)
	middle = Array{AbstractModel, 1}(undef, nBins)
	lower = Array{AbstractModel, 1}(undef, nBins)

	for bin_num in 1:nBins
		# Build models for each layer for this bin
		ignFile = folder_prefix*"/"*data_location*"/bsr_ignmajors_$(bin_num).csv"
		ign, header = readdlm(ignFile, ',', header=true)

		# Read perplex results. Shape (4, 3, n), property, layer, index.
		perplexFile = folder_prefix*"/"*data_location*"/perplex_out_$(bin_num).h5"
		perplexresults = h5read(perplexFile, "results")

		if size(ign,1) != size(perplexresults,3)
			throw(AssertionError("Size of ign does not match size of perplex results from $(fileName)"))
		end

		apply_cracking!(perplexresults, crackProfile, true)

		# Some seismic props may be nan, don't use those.
		okupper = .!isnan.(sum(perplexresults, dims=1)[:,1,:])[:]
		okmiddle = .!isnan.(sum(perplexresults, dims=1)[:,2,:])[:]
		oklower = .!isnan.(sum(perplexresults, dims=1)[:,3,:])[:]
		println("Nan seismic properties found in $(sum(.!okupper)) upper samples, $(sum(.!okmiddle)) middle samples, $(sum(.!oklower)) lower samples.")

		upper[bin_num] = modelType(ign, Array(perplexresults[:,1,:]'))
		middle[bin_num] = modelType(ign, Array(perplexresults[:,2,:]'))
		lower[bin_num] = modelType(ign, Array(perplexresults[:,3,:]'))
	end

	return ModelCollection(Dict(LAYER_NAMES[1]=>upper, LAYER_NAMES[2]=>middle, LAYER_NAMES[3]=>lower), bins, nBins)
end
