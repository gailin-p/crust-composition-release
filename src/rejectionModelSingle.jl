using Distributions
using Random

"""
	SingleVarRejectionModel
Version of rejection model that only uses one seismic var (vp or vs).
TODO: some duplicated code, could combine as subclass of rejection model with some refactoring
"""
abstract type SingleVarRejectionModel <: AbstractModel end

struct VpRejectionModel <: SingleVarRejectionModel
	comp::Array{Float64, 2} # cols according to PERPLEX_ELEMENTS
	seismic::Array{Float64, 2} # cols index, vp
	er::MvNormal
end

struct VsRejectionModel <: SingleVarRejectionModel
	comp::Array{Float64, 2} # cols according to PERPLEX_ELEMENTS
	seismic::Array{Float64, 2} # cols index, vs
	er::MvNormal
end

"""
	VpRejectionModel
Shared constructor takes ign and seismic where cols are index, rho, vp, vpvs
"""
function VpRejectionModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	# error distribution
	dat, header = readdlm("data/adjustedPerplexRes.csv", ',', header=true);
	header = header[:]

	# No eclogites plz
	eclogites = [9,10,11,14,16,19]
	filter = [! (s in eclogites) for s in dat[:,1]]
	dat = dat[filter, :]

	# Interested in non-nan differences between actual and Perplex data
	ers = dat[:,findfirst(isequal("perplex vp"), header)] .- dat[:,findfirst(isequal("dabie vp"), header)]
	ers = ers[.!(isnan.(ers))]
	N_er = fit(MvNormal, ers')

	# Systematic bias may be caused by cracks, deformation caused by transport to surface -- we don't want it.
	zero_mean = MvNormal(zeros(1), N_er.Σ)

	snew = fill(NaN, (size(seismic,1),2))
	snew[:,1] .= seismic[:,1]
	snew[:,2] .= seismic[:,3] # vp

	return VpRejectionModel(ign, snew, zero_mean)
end

function VsRejectionModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	# error distribution
	dat, header = readdlm("data/adjustedPerplexRes.csv", ',', header=true);
	header = header[:]

	# No eclogites plz
	eclogites = [9,10,11,14,16,19]
	filter = [! (s in eclogites) for s in dat[:,1]]
	dat = dat[filter, :]

	# Interested in non-nan differences between actual and Perplex data
	perplex_vs = dat[:,findfirst(isequal("perplex vp"), header)] ./ dat[:,findfirst(isequal("perplex vp/vs"), header)]
	dabie_vs = dat[:,findfirst(isequal("dabie vp"), header)] ./ dat[:,findfirst(isequal("dabie vp/vs"), header)]

	ers = perplex_vs .- dabie_vs
	ers = ers[.!(isnan.(ers))]

	N_er = fit(MvNormal, ers')

	# Systematic bias may be caused by cracks, deformation caused by transport to surface -- we don't want it.
	zero_mean = MvNormal(zeros(1), N_er.Σ)

	snew = fill(NaN, (size(seismic,1),2))
	snew[:,1] .= seismic[:,1]
	snew[:,2] .= seismic[:,3] ./ seismic[:,4] # vs = vp / vpvs

	return VsRejectionModel(ign, snew, zero_mean)
end

function estimateComposition(model::VsRejectionModel,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})
	return estimateComposition(model, vp ./ vpvs)
end

function estimateComposition(model::VpRejectionModel,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})
	return estimateComposition(model, vp)
end

function estimateComposition(model::SingleVarRejectionModel,
	svar::Array{Float64,1})

	toreturn = fill(NaN, (length(svar), size(model.comp,2)))

	for samp in 1:length(svar)
		for i in 1:100000 # give up eventually if really unusual seismic props
			idx = sample(1:size(model.comp,1))
			s = model.seismic[idx,2]
			p = pdf(model.er, [s - svar[samp]])
			maxp = pdf(model.er, model.er.μ) # make alg more efficient by limiting rand vals from 0 to max p
			if rand()*maxp < p
				toreturn[samp,:] .= model.comp[idx,:]
				#println("Selected sample $(model.comp[idx])")
				break
			end
			if i == 100000
				#println("Did not find sample ")
			end
		end
	end

	return toreturn, zeros(size(toreturn))
end

# Just returns one elemnt of interest
function resultSize(model::SingleVarRejectionModel)
	return length(PERPLEX_ELEMENTS) ##
end
