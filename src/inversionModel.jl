"""
Define inversion model type
and functions to estimate composition for given seismic properties for a layer
"""

using MultivariateStats
using StatGeochem
using StatsBase
using Interpolations
using ProgressMeter
using DelimitedFiles
using HDF5
using MAT 

include("../bin.jl") # TODO use Brenhin's version (bin_bsr with nresamples=0)
include("crustDistribution.jl")
include("config.jl")
include("utilities.jl")

prop_labels = ["rho,kg/m3","vp,km/s","vp/vs"] # properties in this order

"""
Information for normalizing samples 
"""
struct Norm 
	std::Float64
	mean::Float64
end 

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

function normalize(norm::Norm, arr::Array)
	return (arr .- norm.mean) ./ norm.std
end

"""
	RangeModel
A model which inverts by returning the mean and std of the compositions of all 
perplex samples within some tolerance of the test sample rho, vp, and vpvs.  
Compare to InversionModel, which uses a binned relationship between PCA and comp
"""
mutable struct RangeModel <: AbstractModel
	# Sorted data. Lookups are sorted arrays, perms are indices back to unsorted data. 
	perms::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}
	lookups::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}
	# Size of bin for each seismic data type 
	errors::Tuple{Float64,Float64,Float64}
	# Original data 
	seismic::Array{Float64,2} # shape (n,4). Columns are index, rho, vp, vpvs 
	comp::Array{Float64,2} # shape (n, n_compositions + 1). Columns are index, wt % major elements 
end 

function resultSize(model::RangeModel)
	return size(model.comp,2)-1
end

"""
	Explicitly set the error ranges for a RangeModel 
"""
function RangeModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2}, errors::Tuple{Float64,Float64,Float64})
	model = RangeModel(ign, seismic)
	model.errors = errors .* (nanmean(model.lookups[1]), nanmean(model.lookups[2]), nanmean(model.lookups[2]))
	return model 
end 

"""
    RangeModel(ign, seismic)
`ign` is size (n, 2) where columns are index, element of interest 
`seismic` is size (n, 4) where columns are index, rho, vp, vp/vs *at this layer!*
"""
function RangeModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	if (seismic[:,1] != 1:size(seismic,1)) | (ign[:,1] != 1:size(seismic,1))
		throw(AssertionError("Expected indicies of seismic and ign data to be 1 through n, with none missing"))
	end 

	# Don't consider samples with any nan values 
	good = ((.~ isnan.(sum(ign,dims=2))) .& (.~(isnan.(sum(seismic,dims=2)))))[:] 
	ign = ign[good,:]
	seismic = seismic[good,:]

	rho_perm = sortperm(seismic[:,2])
	rho_lookup = seismic[:,2][rho_perm]

	vp_perm = sortperm(seismic[:,3])
	vp_lookup = seismic[:,3][vp_perm]

	vpvs_perm = sortperm(seismic[:,4])
	vpvs_lookup = seismic[:,4][vpvs_perm]

	# Default errors TODO I think defaults from igncn1 are too small 
	#all_errors = matread(IGN_FILE)["err2srel"]
	#errors = (all_errors["Rho"]/2, all_errors["Vp"]/2, all_errors["Vp"]/2) # TODO this is not right, but really we want std(vp/vs), and there's no easy way to get that out of std(vp) and std(vs), see https://stats.stackexchange.com/questions/58800/what-is-the-mean-and-standard-deviation-of-the-division-of-two-random-variables
	#errors = errors .* (mean(rho_lookup), mean(vp_lookup), mean(vpvs_lookup))
	
	errors = (.05, .05, .05) .* (nanmean(rho_lookup), nanmean(vp_lookup), nanmean(vpvs_lookup))

	return RangeModel((rho_perm, vp_perm, vpvs_perm), (rho_lookup, vp_lookup, vpvs_lookup), errors, seismic, ign)
end 

"""
	estimateComposition
naiive implementation: search in sorted seismic lists for points within threshold. see: searchsortedfirst
Alternative, if that's too slow: try sorting in 3d space, see GeometricalPredicates 
"""
function estimateComposition(model::RangeModel, 
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})

	if (length(rho) != length(vp)) | (length(vp) != length(vpvs))
		throw(AssertionError("Lengths of seismic datasets don't match"))
	end  

	results = fill(NaN, (length(rho), size(model.comp,2)-1))
	result_stds = fill(NaN, (length(rho), size(model.comp,2)-1))
	#results_indices = fill([], length(rho)) # list of lists of indices for each sample 

	numMatched = RunningMean()
	
	for test_i in 1:length(rho)
		#println("test $test_i")
		# find bottom and top of tolerence
		dat = (rho[test_i], vp[test_i], vpvs[test_i])
		bottoms = dat .- model.errors 
		tops = dat .+ model.errors 

		#println("Gathering samples in seismic range from $bottoms to $tops")

		# Special case: too small or too large, search won't return valid range 
		too_small = [tops[i] < model.lookups[i][1] for i in 1:3] # any top value smaller than bottom 
		too_large = [bottoms[i] > model.lookups[i][end] for i in 1:3] # any bottom value larger than max
		if sum(too_small) + sum(too_large) > 0 
			results[test_i] = NaN 
			result_stds[test_i] = NaN 
			#println("Warning: test sample outside of training data range.")
			continue 
		end 

		# Look up in lookup, then match back to index using perm 
		starts = [searchsortedfirst(model.lookups[i], bottoms[i]) for i in 1:3]
		ends = [searchsortedfirst(model.lookups[i], tops[i]) for i in 1:3]
		ends = ends .- 1 # searchsorted returns index of first value *larger than* query 
		# Match back to indices 
		indices = [Set(model.perms[i][starts[i]:ends[i]]) for i in 1:3]

		# find indices agreed upon by all types of seismic data 
		# TODO make this faster if we end up using big bins, for now assume n^2 is fine 
		agreed = Array{Int, 1}(undef, 0)
		for i in indices[1]
			if (i in indices[2]) & (i in indices[3])
				push!(agreed, Int(i))
			end 
		end 
		#println("agree on $agreed")
		# How many were matched? 
		mean!(numMatched, length(agreed))

		# look up compositions of agreed upon indices 
		compositions = model.comp[agreed, 2:end] # First is just index 
		#results_indices[test_i] = model.comp[agreed, 1]
		results[test_i,:] = mean(compositions, dims=1) 
		result_stds[test_i,:] = std(compositions, dims=1)


		# Sanity check: mean seismic data for agreed-upon indices better be within error of test val.
		matched_rho = nanmean(model.seismic[agreed, 2])
		matched_vp = nanmean(model.seismic[agreed, 3])
		matched_vpvs = nanmean(model.seismic[agreed, 4])
		#println("matched rho, vp, vpvs = $matched_rho, $matched_vp, $matched_vpvs")
		if (matched_rho > tops[1]) | (matched_rho < bottoms[1]) | 
			(matched_vp > tops[2]) | (matched_vp < bottoms[2]) | 
			(matched_vpvs > tops[3]) | (matched_vpvs < bottoms[3]) 
			throw(AssertionError("Model matched samples which don't have seismic data close to test."))
		end 
	end

	# println("reporting from estimation function: 
	# 	Mean matched samples $(numMatched.m) for $(numMatched.n) samples tried.
	# 	$(sum(.! isnan.(results))) not NaN results out of $(length(rho)) inputs")
	if numMatched.m < 10
		println("Warning: matched on average less than 10 perplex samples per test sample.")
	end 

	return results, result_stds #, results_indices

end

# Just returns one elemnt of interest 
function resultSize(model::InversionModel)
	return 1
end

function resultSize(model::BadModel)
	return 1 
end 

"""
    InversionModel(ign, seismic)
`ign` is size (n, 2) where columns are index, element of interest 
`seismic` is size (n, 4) where columns are index, rho, vp, vp/vs *at this layer!*
"""
function InversionModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	# Check that indices align 
	if ign[:,1] != seismic[:,1]
		#println(hcat(ign[1:100,1], seismic[1:100,1]))
		throw(ArgumentError("Indices in seismic data do not match those in composition data"))
	end

	# Only use SiO2 
	ign = ign[:,1:2]

	# Normalize data 
	rhoNorm = Norm(nanstd(seismic[:,2]), nanmean(seismic[:,2]))
	vpNorm = Norm(nanstd(seismic[:,3]), nanmean(seismic[:,3]))
	vpvsNorm = Norm(nanstd(seismic[:,4]), nanmean(seismic[:,4]))
    samples = vcat(normalize(rhoNorm,seismic[:,2])', 
    			   normalize(vpNorm,seismic[:,3])',
    			   normalize(vpvsNorm,seismic[:,4])')

    # Discard NaN values 
    test = .!(isnan.(samples[1,:]))
    samples = samples[:,test] 

    # Do we have enough valid perplex results in this bin? 
    if size(samples,2) < 5 # Todo is this a reasonable cutoff? 
    	println("Warning: not enough samples to build model for this layer and geotherm.")
    	return BadModel()
    end

    # PCA
    pca = fit(PCA, samples)
    pca_samples = transform(pca, samples)[1,:]

    # Bin  
    # Look at your independent (seismic) variable of choice
    xmin = percentile(pca_samples, .5)
    xmax = percentile(pca_samples,99.5)

    # Bin your element of interest as a function of your seismic variable
    # center (seismic), mean (element of interest eg SiO2), error
    # TODO not using oversampling ratio right now. what's correct? 
    c, m, e = bin(pca_samples, ign[:,2][test], xmin, xmax, 40)

    itp = interpolate(m, BSpline(Linear())) # Linear interpolation of bin means
    itp = scale(itp, c) # Scale by bin centers
    itp = extrapolate(itp, Flat()) # outside bin means, just give the upper or lower bin mean. TODO gailin is this the right choice?

    # Interpolate errors
    itp_e = interpolate(e, BSpline(Linear())) # Linear interpolation of bin means
    itp_e = scale(itp_e, c) # Scale by bin centers
    itp_e = extrapolate(itp_e, Flat()) # TODO gailin - why is it ok to interpolate uncertainty?

    return InversionModel(rhoNorm, vpNorm, vpvsNorm, pca, itp, itp_e, c, m, e, pca_samples, ign[:,2][test])
end

"""
    estimateComposition(model, rho, vp, vpvs)

Return an array of the same length as the input params 
"""
function estimateComposition(model::InversionModel, 
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})
	if (length(rho) != length(vp)) | (length(vp) != length(vpvs))
		throw(AssertionError("Lengths of seismic datasets don't match"))
	end

	# Normalize and run PCA on input data 
    samples = vcat(normalize(model.rhoNorm,rho)', 
    			   normalize(model.vpNorm,vp)',
    			   normalize(model.vpvsNorm,vpvs)')

    # Discard NaN values 
    test = .!(isnan.(samples[1,:]))
    samples = samples[:,test] 

    # PCA
    pca_samples = transform(model.pca, samples)[1,:]

    # Run on test data 
    m_itp = model.interp(pca_samples)

    # Interpolate errors
    e_itp = model.error_interp(pca_samples)

    return (m_itp, e_itp)
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

	# Input and result data for all geotherm bins 
	nresults = resultSize(models.models["upper"][1])
	results = fill(0.0, (length(rho), nresults))
	errors = fill(0.0, (length(rho), nresults))
	done = fill(0, length(rho))

	@showprogress "Running samples" for (i, bin_bottom) in enumerate(models.bins[1:end-1])
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
		results[test,:] = means
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
function makeModels(data_location::String; resamplePerplex::Bool=false, modelType::DataType=InversionModel)
	# Elements in ign files: index, elts, geotherm, layers. TODO add header to CSV created by resampleEarthChem
	elements = PERPLEX_ELEMENTS

	# Build models for every (geotherm bin)/layer combo 
	fileNames = filter(x->contains(x,"perplex_out_"), readdir("data/$data_location"))
	nBins = length(fileNames)
	if nBins == 0 
		throw(AssertionError("No perplex results found for data directory $data_location"))
	end 
	bins = crustDistribution.binBoundaries(nBins)

	# Models for each layer. upper[i] is model for i-th geotherm bin. 
	upper = Array{AbstractModel, 1}(undef, nBins)
	middle = Array{AbstractModel, 1}(undef, nBins)
	lower = Array{AbstractModel, 1}(undef, nBins)

	@showprogress 1 "Building models" for bin_num in 1:nBins
		# Build models for each layer for this bin 
		ignFile = "data/"*data_location*"/bsr_ignmajors_$(bin_num).csv"
		ign = readdlm(ignFile, ',')

		# Read perplex results. Shape (4, 3, n), property, layer, index. 
		perplexFile = "data/"*data_location*"/perplex_out_$(bin_num).h5"
		perplexresults = h5read(perplexFile, "results")

		# Resample perplex? 
		if resamplePerplex
			throw(AssertionError("resampling not implemented yet srry"))
			bsrresample() # data, sigma, nrows 
		end 

		if size(ign,1) != size(perplexresults,3)
			throw(AssertionError("Size of ign does not match size of perplex results from $(fileName)"))
		end 

		# Build models. use SiO2 (index 2 in ign). 
		upper[bin_num] = modelType(ign, Array(perplexresults[:,1,:]'))
		middle[bin_num] = modelType(ign, Array(perplexresults[:,2,:]'))
		lower[bin_num] = modelType(ign, Array(perplexresults[:,3,:]'))
	end 

	return ModelCollection(Dict(LAYER_NAMES[1]=>upper, LAYER_NAMES[2]=>middle, LAYER_NAMES[3]=>lower), bins, nBins)
end 

