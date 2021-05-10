using Distributions

"""
	RangeModel
A model which inverts by returning the mean and std of the compositions of all 
perplex samples within some tolerance of the test sample rho, vp, and vpvs.  
Compare to InversionModel, which uses a binned relationship between PCA and comp
"""
struct RangeModel <: AbstractModel
	# Sorted data. Lookups are sorted arrays, perms are indices back to unsorted data. 
	# perms::Tuple{Array{Int,1}, Array{Int,1}, Array{Int,1}}
	# lookups::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}} # rho, vp, vpvs 
	perms::Array{Int,2} # inversion var, perm 
	lookups::Array{Float64,2} # inversion var, sample. this order means all of each var are together in memory
	# Size of bin for each inversion var
	errors::Array{Float64}
	# Original data 
	seismic::Array{Float64,2} # shape (n,4). Columns are index, rho, vp, vpvs 
	comp::Array{Float64,2} # shape (n, n_compositions + 1). Columns are index, wt % major elements 
	#mean::Bool # Take mean of composition samples matching each test sample? 
	#use::Tuple{Bool, Bool, Bool} # Which of vp, rho, vpvs to use? 
	er::MvNormal
end 

"""
	A version of range model using only vp 
	Assumes input `seismic` has cols index, rho, vp, vpvs 

"""
function VpRangeModel(ign, seismic)
	# error distribution 
	dat, header = readdlm("data/adjustedPerplexRes.csv", ',', header=true);
	header = header[:]
	# Interested in non-nan differences between actual and Perplex data 
	ers = fill(NaN, (size(dat)[1],1));
	ers[:,1] .= dat[:,findfirst(isequal("perplex vp"), header)] .- dat[:,findfirst(isequal("dabie vp"), header)]
	ers = ers[.!(isnan.(sum(ers, dims=2)))[:],:]; 
	N_er = fit(MvNormal, ers')

	# Bias may be caused by cracks, deformation caused by transport to surface -- we don't want it. 
	zero_mean = MvNormal(zeros(1), N_er.Σ)

	return RangeModel(ign, hcat(seismic[:,1], seismic[:,3]))
end

"""
	A version of range model using only vp and rho 
	Assumes input `seismic` has cols index, rho, vp, vpvs 
"""
function VpVsRhoRangeModel(ign, seismic)
	return RangeModel(ign, hcat(seismic[:,1], seismic[:,3]))
end


"""
	RangeModel(ign, seismic)
Default constructor for a rho, vp, vpvs model. 
"""
function RangeModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2}; relerr=.01)
	# error distribution 
	dat, header = readdlm("data/adjustedPerplexRes.csv", ',', header=true);
	header = header[:]
	# Interested in non-nan differences between actual and Perplex data 
	ers = fill(NaN, (size(dat)[1],3));
	ers[:,2] .= dat[:,findfirst(isequal("perplex vp"), header)] .- dat[:,findfirst(isequal("dabie vp"), header)]
	ers[:,3] .= dat[:,findfirst(isequal("perplex vp/vs"), header)] .- dat[:,findfirst(isequal("dabie vp/vs"), header)]
	ers[:,1] .= dat[:,findfirst(isequal("perplex rho"), header)] .- (1000 .* dat[:,findfirst(isequal("dabie rho"), header)]);
	ers = ers[.!(isnan.(sum(ers, dims=2)))[:],:]; 
	N_er = fit(MvNormal, ers')

	# Bias may be caused by cracks, deformation caused by transport to surface -- we don't want it. 
	zero_mean = MvNormal(zeros(3), N_er.Σ)

	return RangeModel(ign, seismic, zero_mean, relerr=relerr)
end 

"""
    RangeModel(ign, seismic, uncertainty)

Shared constructor for any range model. 
`ign` is size (n, x) where columns are cols in bsr_ignmajors_[bin].csv 
`seismic` is size (n, 4) where columns are index, inversion vars *at this layer!*
`uncertainty` is a distribution, ie, I can call rand(uncertainty) and get a sample of perplex seismic error 
"""
function RangeModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2}, uncertanty::Distribution; relerr=.01, uncertain=true)
	if (seismic[:,1] != 1:size(seismic,1)) | (ign[:,1] != 1:size(seismic,1))
		throw(AssertionError("Expected indicies of seismic and ign data to be 1 through n, with none missing"))
	end 

	# Don't consider samples with any nan values in major elts. 
	good = ((.~ isnan.(sum(ign[:,2:1+length(COMPOSITION_ELEMENTS)],dims=2))) # Major elts ok 
		.& (.~(isnan.(sum(seismic,dims=2)))))[:] # Seismic ok 
	println("Using $(sum(good)) compositions in inversion.")
	ign = ign[good,:]
	seismic = seismic[good,:]

	nvars = size(seismic,2) - 1 # num columns minus index column 
	nsamples = size(seismic,1)

	perms = fill(1, (nvars, nsamples))
	lookups = fill(NaN, (nvars, nsamples))

	for i in 1:nvars 
		perm = sortperm(seismic[:,i+1])
		perms[i,:] .= perm 
		lookups[i,:] .= seismic[perm,i+1]
	end

	search_range = relerr .* nanmean(lookups, dims=2)

	return RangeModel(perms, lookups, search_range, seismic, ign, uncertanty)
end 

function estimateComposition(model::RangeModel, dat::Array{Float64,1})
	newdat = Array{Float64,2}(undef,(length(dat),1))
	newdat[:,1] .= dat
	return estimateComposition(model, newdat)
end

function estimateComposition(model::RangeModel, 
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1}) 

	if (length(rho) != length(vp)) | (length(vp) != length(vpvs))
		throw(AssertionError("Lengths of seismic datasets don't match"))
	end 

	estimateComposition(model, hcat(rho, vp, vpvs))
end 

"""
	estimateComposition
naiive implementation: search in sorted seismic lists for points within threshold. see: searchsortedfirst
Alternative, if that's too slow: try sorting in 3d space, see GeometricalPredicates 

inversionDat::Array{Float64,2} Rows = samples to invert, cols = inversion vars. 
	Number of inversion vars must equal number of inversion vars in model. 


Proof: works as long as the min distance inverted sample is within relerr % of input sample. 
TODO: make a heatmap of what parts of vp/rho space fufil this assumption for a given perplex dset 
"""
function estimateComposition(model::RangeModel, inversionDat::Array{Float64,2}) 

	if size(inversionDat,2) != size(model.perms,1)
		throw(AssertionError("Number of inversion vars don't match"))
	end  
	nvars = size(inversionDat,2)

	results = fill(NaN, (size(inversionDat,1), size(model.comp,2)))
	result_stds = fill(NaN, (size(inversionDat,1), size(model.comp,2)))
	#results_indices = fill([], length(rho)) # list of lists of indices for each sample 

	numMatched = RunningMean()
	
	for test_i in 1:size(inversionDat,1)
		for attempt in 1:50
			#println("test $test_i")
			# find bottom and top of tolerence
			dat = inversionDat[test_i,:] .+ rand(model.er)
			bottoms = dat .- model.errors 
			tops = dat .+ model.errors 

			#println("Gathering samples in seismic range from $bottoms to $tops")

			# Special case: too small or too large, search won't return valid range 
			too_small = tops .< model.lookups[:,1] # any top value smaller than bottom 
			too_large = bottoms .> model.lookups[:,end] # any bottom value larger than max
			if sum(too_small) + sum(too_large) > 0 
				results[test_i] = NaN 
				result_stds[test_i] = NaN 
				#println("Warning: test sample outside of training data range.")
				continue 
			end 

			# Look up in lookup, then match back to index using perm 
			starts = [searchsortedfirst(model.lookups[i,:], bottoms[i]) for i in 1:nvars]
			ends = [searchsortedlast(model.lookups[i,:], tops[i]) for i in 1:nvars]
			# Match back to indices 
			# lengths = (ends .- starts) .+ 1
			# if minimum(lengths) == 0 
			# 	results[test_i,:] .= NaN
			#  	result_stds[test_i,:] .= NaN
			#  	continue 
			# end

			indices = [Set(model.perms[i,starts[i]:ends[i]]) for i in 1:nvars]

			#agreed = Array{Int, 1}(undef, 0) ## indices of comps within error % of seismic sample 
			#distances = Array{Float64, 1}(undef, 0) ## L1 norm distance of each agreed comp to seismic sample 
			idx = 0
			best = Inf 
			for i in indices[1]
				if nvars > 1 
					check = all([i in indices[j] for j in 2:nvars])
					if check
					#push!(agreed, i)
						distance = abs(sum((model.seismic[i,2:end] .- dat) ./ dat)) # L1 norm of percent difference 
						if distance < best 
							idx = i 
							best = distance 
						end
					end
				else 
					distance = abs(sum((model.seismic[i,2:end] .- dat) ./ dat)) # L1 norm of percent difference 
					if distance < best 
						idx = i 
						best = distance 
					end
				end 
			end 

			if idx == 0 
				results[test_i,:] .= NaN
			 	result_stds[test_i,:] .= NaN
			 	continue 
			else 
				results[test_i,:] = model.comp[idx, :]
				result_stds[test_i,:] .= -1
				break
			end

			# use = argmin(lengths) # use var with shortest matched list for SPEEDY 

			# # find index of closest 
			# indices = model.perms[use,starts[use]:ends[use]]
			# best = Inf
			# idx = 0
			# for possible in indices
			# 	distance = sum((model.seismic[possible,2:end] .- dat) ./ dat) # L1 norm of percent difference
			# 	if abs(distance) < best 
			# 		best = abs(distance)
			# 		idx = possible
			# 	end
			# end 

			
		end

	end

	return results, result_stds #, results_indices

end


function resultSize(model::RangeModel)
	return size(model.comp,2)
end

