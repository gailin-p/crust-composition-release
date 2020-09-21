"""
Implement InversionModel interface for a version of RangeModel which only uses Vp (not Vp/Vs or Rho)
"""

include("inversionModel.jl")


"""
	VpModel
A model which inverts by returning the mean and std of the compositions of all 
perplex samples within some tolerance of the test sample rho, vp, and vpvs.  
Compare to InversionModel, which uses a binned relationship between PCA and comp
"""
mutable struct VpModel <: AbstractModel
	# Sorted data. Lookups is sorted arrays, perms is indices back to unsorted data. 
	perm::Array{Float64,1}
	lookup::Array{Float64,1} 
	# Size of bin as Vp 
	bin::Float64
	# Original data 
	seismic::Array{Float64,2} # shape (n,4). Columns are index, rho, vp, vpvs 
	comp::Array{Float64,2} # shape (n, n_compositions + 1). Columns are index, wt % major elements 
end 

function resultSize(model::VpModel)
	return size(model.comp,2)-1
end

"""
	Explicitly set the error ranges for a VpModel 
Only use errors[2], the vp error. 
"""
function setError(model::VpModel, errors::Tuple{Float64,Float64,Float64})
	model.bin = errors[2] * nanmean(model.lookup)
	return model 
end 


"""
    VpModel(ign, seismic)
`ign` is size (n, 2) where columns are index, element of interest 
`seismic` is size (n, 4) where columns are index, rho, vp, vp/vs *at this layer!* 
All of rho, vp/vs will be ignored. 
"""
function VpModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	if (seismic[:,1] != 1:size(seismic,1)) | (ign[:,1] != 1:size(seismic,1))
		throw(AssertionError("Expected indicies of seismic and ign data to be 1 through n, with none missing"))
	end 

	# Don't consider samples with any nan values 
	good = ((.~ isnan.(sum(ign,dims=2))) .& (.~(isnan.(sum(seismic,dims=2)))))[:] 
	ign = ign[good,:]
	seismic = seismic[good,:]

	vp_perm = sortperm(seismic[:,3])
	vp_lookup = seismic[:,3][vp_perm]

	# Default bin sizes 
	bin = .05 * nanmean(vp_lookup)

	return VpModel(vp_perm, vp_lookup, bin, seismic, ign)
end 

"""
	estimateComposition
Only vp used. 
"""
function estimateComposition(model::VpModel, 
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})

	results = fill(NaN, (length(vp), size(model.comp,2)-1))
	result_stds = fill(NaN, (length(vp), size(model.comp,2)-1))

	numMatched = RunningMean()
	
	for test_i in 1:length(vp)
		#println("test $test_i")
		# find bottom and top of tolerence
		sample = vp[test_i]
		bottom = sample - model.bin
		top = sample + model.bin

		# Special case: too small or too large, search won't return valid range 
		if (top < model.lookup[1]) | (bottom > model.lookup[end])
			results[test_i] = NaN 
			result_stds[test_i] = NaN 
			#println("Warning: test sample outside of training data range.")
			continue 
		end 

		# Look up in lookup, then match back to index using perm 
		start = searchsortedfirst(model.lookup, bottom) 
		ends = searchsortedfirst(model.lookup, top) 
		ends = ends .- 1 # searchsorted returns index of first value *larger than* query 
		# Match back to indices 
		agreed = Int.(Array(model.perm[start:ends]))

		#println("agree on $agreed")
		# How many were matched? 
		mean!(numMatched, length(agreed))

		# look up compositions of agreed upon indices 
		compositions = model.comp[agreed, 2:end] # First is just index 
		
		if length(agreed)==0
			results[test_i,:] .= NaN
			result_stds[test_i,:] .= NaN
		else
			# Select one of matching samples 
			idx = rand(1:length(agreed))
			results[test_i,:] = compositions[idx, :]
			result_stds[test_i,:] .= -1
			
			# Alternative: Mean sample 
			#results[test_i,:] = mean(compositions, dims=1) 
			#result_stds[test_i,:] = std(compositions, dims=1)
		end

	end

	if numMatched.m < 10
		println("Warning: matched on average less than 10 perplex samples per test sample.")
	end 

	return results, result_stds 
end 

