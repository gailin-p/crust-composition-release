"""
Implement InversionModel interface for a version of RangeModel which only uses Vp and Rho (not Vp/Vs)
"""

include("inversionModel.jl")


"""
	VpRhoModel
A model which inverts by returning the mean and std of the compositions of all 
perplex samples within some tolerance of the test sample rho, vp, and vpvs.  
Compare to InversionModel, which uses a binned relationship between PCA and comp
"""
mutable struct VpRhoModel <: AbstractModel
	# Sorted data. Lookups are sorted arrays, perms are indices back to unsorted data. 
	perms::Tuple{Array{Float64,1}, Array{Float64,1}}
	lookups::Tuple{Array{Float64,1}, Array{Float64,1}} # rho, vp
	# Size of bin for each seismic data type 
	errors::Tuple{Float64,Float64}
	# Original data 
	seismic::Array{Float64,2} # shape (n,4). Columns are index, rho, vp 
	comp::Array{Float64,2} # shape (n, n_compositions + 1). Columns are index, wt % major elements 
	mean::Bool # use mean instead of random of matching composition samples 
end 

function resultSize(model::VpRhoModel)
	return size(model.comp,2)-1
end

"""
	Explicitly set the error ranges for a VpRhoModel 
Only use errors[2], the vp error. 
"""
function setError(model::VpRhoModel, errors::Tuple{Float64,Float64,Float64})
	model.errors = (errors[1]*nanmean(model.lookups[1]), errors[2]*nanmean(model.lookups[2]))
	return model 
end 

"""
	Explicitly set the error ranges for a RangeModel 
"""
function setMean(model::VpRhoModel, usemean::Bool)
	model.mean = usemean
	return model 
end 


"""
    VpRhoModel(ign, seismic)
`ign` is size (n, 2) where columns are index, element of interest 
`seismic` is size (n, 4) where columns are index, rho, vp, vp/vs *at this layer!* 
All of rho, vp/vs will be ignored. 
"""
function VpRhoModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	if (seismic[:,1] != 1:size(seismic,1)) | (ign[:,1] != 1:size(seismic,1))
		throw(AssertionError("Expected indicies of seismic and ign data to be 1 through n, with none missing"))
	end 

	# Don't consider samples with any nan values 
	good = ((.~ isnan.(sum(ign,dims=2))) .& (.~(isnan.(sum(seismic,dims=2)))))[:] 
	ign = ign[good,:]
	seismic = seismic[good,:]

	vp_perm = sortperm(seismic[:,3])
	vp_lookup = seismic[:,3][vp_perm]

	rho_perm = sortperm(seismic[:,2])
	rho_lookup = seismic[:,2][rho_perm]

	# Default bin sizes 
	bin = (.05 * nanmean(vp_lookup), .05 * nanmean(rho_lookup))

	return VpRhoModel((rho_perm, vp_perm), (rho_lookup, vp_lookup), bin, seismic, ign, false)
end 

"""
	estimateComposition
Only vp used. 
"""
function estimateComposition(model::VpRhoModel, 
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1}) 

	if length(rho) != length(vp)
		throw(AssertionError("Lengths of seismic datasets don't match"))
	end  

	results = fill(NaN, (length(rho), size(model.comp,2)-1))
	result_stds = fill(NaN, (length(rho), size(model.comp,2)-1))
	#results_indices = fill([], length(rho)) # list of lists of indices for each sample 

	numMatched = RunningMean()
	
	for test_i in 1:length(rho)
		#println("test $test_i")
		# find bottom and top of tolerence
		dat = (rho[test_i], vp[test_i])
		bottoms = dat .- model.errors 
		tops = dat .+ model.errors 

		#println("Gathering samples in seismic range from $bottoms to $tops")

		# Special case: too small or too large, search won't return valid range 
		too_small = [tops[i] < model.lookups[i][1] for i in 1:2] # any top value smaller than bottom 
		too_large = [bottoms[i] > model.lookups[i][end] for i in 1:2] # any bottom value larger than max
		if sum(too_small) + sum(too_large) > 0 
			results[test_i] = NaN 
			result_stds[test_i] = NaN 
			continue 
		end 

		# Look up in lookup, then match back to index using perm 
		starts = [searchsortedfirst(model.lookups[i], bottoms[i]) for i in 1:2]
		ends = [searchsortedfirst(model.lookups[i], tops[i]) for i in 1:2]
		ends = ends .- 1 # searchsorted returns index of first value *larger than* query 
		# Match back to indices 
		indices = [Set(model.perms[i][starts[i]:ends[i]]) for i in 1:2]

		agreed = Array{Int, 1}(undef, 0)
		for i in indices[1]
			if (i in indices[2])
				push!(agreed, Int(i))
			end 
		end 

		mean!(numMatched, length(agreed))

		# look up compositions of agreed upon indices 
		compositions = model.comp[agreed, 2:end] # First is just index 
		
		if length(agreed)==0
			results[test_i,:] .= NaN
			result_stds[test_i,:] .= NaN
		else 
			if !model.mean
				# Select one of matching samples
				idx = rand(1:length(agreed))
				results[test_i,:] = compositions[idx, :]
				result_stds[test_i,:] .= -1
			else 
				# Alternative: Mean sample 
				results[test_i,:] = mean(compositions, dims=1) 
				result_stds[test_i,:] = std(compositions, dims=1)
			end 
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

