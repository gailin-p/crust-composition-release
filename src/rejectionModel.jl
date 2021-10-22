using Distributions
using Random

"""
	RejectionModel
A model which inverts by using rejection sampling to sample from posterior
distribution P(comp | seismic) for seismic prop at every location
"""
struct RejectionModel <: AbstractModel
	comp::Array{Float64, 2} # cols according to PERPLEX_ELEMENTS
	seismic::Array{Float64, 2} # cols index, rho, vp, vpvs
	er::MvNormal
end

function RejectionModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	# error distribution
	dat, header = readdlm("data/adjustedPerplexRes.csv", ',', header=true);
	header = header[:]

	# No eclogites plz
	eclogites = [9,10,11,14,16,19]
	filter = [! (s in eclogites) for s in dat[:,1]]
	dat = dat[filter, :]

	# Interested in non-nan differences between actual and Perplex data
	ers = fill(NaN, (size(dat)[1],3));
	ers[:,2] .= dat[:,findfirst(isequal("perplex vp"), header)] .- dat[:,findfirst(isequal("dabie vp"), header)]
	ers[:,3] .= dat[:,findfirst(isequal("perplex vp/vs"), header)] .- dat[:,findfirst(isequal("dabie vp/vs"), header)]
	ers[:,1] .= dat[:,findfirst(isequal("perplex rho"), header)] .- (1000 .* dat[:,findfirst(isequal("dabie rho"), header)]);
	ers = ers[.!(isnan.(sum(ers, dims=2)))[:],:];
	N_er = fit(MvNormal, ers')

	# Systematic bias may be caused by cracks, deformation caused by transport to surface -- we don't want it.
	zero_mean = MvNormal(zeros(3), N_er.Σ) # Divide sigma by 2 for narrower error -- but this does *not* change systematic bias of result. 

	return RejectionModel(ign, seismic, zero_mean)
end

function estimateComposition(model::RejectionModel,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})

	toreturn = fill(NaN, (length(rho), size(model.comp,2)))

	for samp in 1:length(rho)
		for i in 1:100000 # give up eventually if really unusual seismic props
			idx = sample(1:size(model.comp,1))
			s = model.seismic[idx,2:end]
			p = pdf(model.er, s .- [rho[samp],vp[samp],vpvs[samp]])
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
function resultSize(model::RejectionModel)
	return length(PERPLEX_ELEMENTS) ##
end
