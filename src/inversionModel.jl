"""
Define inversion model type
and functions to estimate composition for given seismic properties for a layer
"""

using MultivariateStats
using StatGeochem
using StatsBase
using Interpolations

include("../bin.jl") # TODO use Brenhin's version (bin_bsr with nresamples=0)


prop_labels = ["rho,kg/m3","vp,km/s","vp/vs"] # properties in this order

"""
Information for normalizing samples 
"""
struct Norm 
	std::Float64
	mean::Float64
end 

"""
    InversionModel 
Store information for inverting observed data according to perplex data for samples 
# Includes: PCA, mean/std for normalizing inputs, relationship between PCA and comp
# Different models for each inversion (upper != middle != lower)
"""
struct InversionModel
	# PCA of seismic vars 
	rhoNorm::Norm
	vpNorm::Norm
	vpvsNorm::Norm
	pca::PCA 
	# Binned relationship between PCA and composition
	c # center of bin (seismic PCA)
	m::Array{Float64,1} # mean of bin (elt of interest)
	e::Array{Float64,1} # error of bin (elt of interest)
end

function normalize(norm::Norm, arr::Array)
	return (arr .- norm.mean) ./ norm.std
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

    # PCA
    pca = fit(PCA, samples)
    pca_samples = transform(pca, samples)[1,:]

    # Bin  TODO below copy-pasted
    # Look at your independent (seismic) variable of choice
    xmin = percentile(pca_samples, .5)
    xmax = percentile(pca_samples,99.5)

    # Bin your element of interest as a function of your seismic variable
    # center (seismic), mean (element of interest eg SiO2), error
    # TODO not using oversampling ratio right now. what's correct? 
    c, m, e = bin(pca_samples, ign[:,2][test], xmin, xmax, 40)

    return InversionModel(rhoNorm, vpNorm, vpvsNorm, pca, c, m, e)
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

    # TODO move the interpolation defn into model creation 
    # Interpolate (element of interest) values for each (obs seismic) sample in the dataset
    itp = interpolate(model.m, BSpline(Linear())) # Linear interpolation of bin means
    sitp = scale(itp, model.c) # Scale by bin centers
    esitp = extrapolate(sitp, Flat()) # outside bin means, just give the upper or lower bin mean. TODO gailin is this the right choice?
    # Run on test data 
    m_itp = esitp(pca_samples)

    # Interpolate errors
    itp_e = interpolate(model.e, BSpline(Linear())) # Linear interpolation of bin means
    sitp_e = scale(itp_e, model.c) # Scale by bin centers
    esitp_e = extrapolate(sitp_e, Flat()) # TODO gailin - why is it ok to interpolate uncertainty?
    
    e_itp = esitp_e(pca_samples)

    return (m_itp, e_itp)
end


