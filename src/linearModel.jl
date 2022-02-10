"""
The simplest inversion model I can think of. Fit Vx = c for each composition elt c,
where V is normalized seismic params plus an index row. Then to invert a sample v,
just return the fit for each elt
"""

include("inversionModel.jl")

"""
    LinearModel 
"""
struct  LinearModel <: AbstractModel
	# PCA of seismic vars
	rhoNorm::Norm
	vpNorm::Norm
	vpvsNorm::Norm
	# Fits
	fit::Array{Float64,2} # Dimensions composition types, then fit

	# Original data
	seismic::Array{Float64,2} # Seismic data of perplex samples post-pca. length n
	comp::Array{Float64,2} # composition of element of interest. length n.
end

function LinearModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	# Build A

	# Normalize data
	rhoNorm = Norm(nanstd(seismic[:,2]), nanmean(seismic[:,2]))
	vpNorm = Norm(nanstd(seismic[:,3]), nanmean(seismic[:,3]))
	vpvsNorm = Norm(nanstd(seismic[:,4]), nanmean(seismic[:,4]))
    A = hcat(normalize(rhoNorm,seismic[:,2]),
    			   normalize(vpNorm,seismic[:,3]),
    			   normalize(vpvsNorm,seismic[:,4]),
    			   fill(1,size(seismic)[1]))

    # Discard NaN values
    test = .!(isnan.(sum(A,dims=2)))[:]
    A = A[test,:]
    ign = ign[test,:]

	# Fit lstsq for each element c
	ncomp = size(ign)[2]
	fits = zeros((ncomp, 4))
	for c in 1:ncomp
		ctest = .!(isnan.(ign[:,c]))[:]
		fits[c,:] .= A[ctest,:] \ ign[ctest,c]
	end

	return LinearModel(rhoNorm, vpNorm, vpvsNorm, fits, seismic, ign)

end

function estimateComposition(model::LinearModel,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})

	rhoNorm = normalize(model.rhoNorm,rho)
    vpNorm = normalize(model.vpNorm,vp)
    vpvsNorm = normalize(model.vpvsNorm,vpvs)

    result = fill(NaN, (length(rho), size(model.fit)[1]))
    errs = fill(NaN, (length(rho), size(model.fit)[1]))

    for c in 1:size(model.fit)[1] # how many dependent vars have fits in the model
    	result[:,c] .= (rhoNorm .* model.fit[c,1]) .+
    					(vpNorm .* model.fit[c,2]) .+
    					(vpvsNorm .* model.fit[c,3]) .+ model.fit[c,4]
    end

    result[:,1] .= 0 # first row is sample index, which doesn't mean anything

    return result, errs
end

# Just returns one elemnt of interest
function resultSize(model::LinearModel)
	return length(PERPLEX_ELEMENTS) ## Just to match inverison model
end
