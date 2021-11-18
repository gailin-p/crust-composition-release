using Distributions
using Random
using GaussianMixtures
using MultivariateStats

"""
	RejectionPriorModel
A model which inverts by using rejection sampling to sample from posterior
distribution P(comp | seismic) for seismic prop at every location

Allows selection of random priors with different, randomly selected mixtures of mafic/felsic rocks
where mafic/felsic peaks are found by fitting spatially resampled earthchem prior (RejectionPriorModel.basePrior)

Intended use: sensitivity testing of prior, where newPrior(model) is called
"""
mutable struct RejectionPriorModel <: AbstractModel
	comp::Array{Float64, 2} # cols according to PERPLEX_ELEMENTS
	seismic::Array{Float64, 2} # cols index, rho, vp, vpvs
	er::MvNormal
	basePrior::MixtureModel # Fitted mixture model of two 10-d gaussians fitted to
	currentPrior::MixtureModel # mixture model currently used.
	means::Array{Float64, 1} # Record center of used distributions
end

function newPrior!(models::ModelCollection)
	m = models.models[UPPER][1]
	a = rand()
    currentPrior = MixtureModel(m.basePrior.components, [a, 1-a])
	for layer in values(models.models)
		for m in layer
			newPrior!(m, currentPrior)
		end
	end
end

"""
	Use a passed in prior so change all models in a model collection at once
"""
function newPrior!(model::RejectionPriorModel, newPrior::MixtureModel)
    model.currentPrior = newPrior
	a = newPrior.prior.p[1]
	append!(model.means,
		a*model.currentPrior.components[1].μ[1] + (1-a)*model.currentPrior.components[2].μ[1])
	@info "Prior mean" last(model.means)
end

"""
	Assign random new prior (for testing)
"""
function newPrior!(model::RejectionPriorModel)
	a = rand()
    currentPrior = MixtureModel(model.basePrior.components, [a, 1-a])
	append!(model.means,
		a*model.currentPrior.components[1].μ[1] + (1-a)*model.currentPrior.components[2].μ[1])
	@info "Prior mean" last(model.means)
end

function RejectionPriorModel(ign::Array{Float64, 2}, seismic::Array{Float64, 2})
	# error distribution
	if isfile("data/adjustedPerplexRes.csv")
		dat, header = readdlm("data/adjustedPerplexRes.csv", ',', header=true);
	else
		dat, header = readdlm("../data/adjustedPerplexRes.csv", ',', header=true);
	end
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

	## Prior distribution
	if isfile("resources/bsr_ignmajors_1.csv")
		prior, h = readdlm("resources/bsr_ignmajors_1.csv", ',', header=true);
	else
		prior, h = readdlm("../resources/bsr_ignmajors_1.csv", ',', header=true);
	end
	gm_model = MixtureModel(GMM(2, prior[:,2:11]))
	model = RejectionPriorModel(ign, seismic, zero_mean, gm_model, gm_model, zeros(0))
	# User is responsible for making random prior using newPrior!()

	return model
end

function estimateComposition(model::RejectionPriorModel,
	rho::Array{Float64,1}, vp::Array{Float64,1}, vpvs::Array{Float64,1})

	toreturn = fill(NaN, (length(rho), size(model.comp,2)))

	for samp in 1:length(rho)
		for i in 1:200000 # give up eventually if really unusual seismic props
			idx = sample(1:size(model.comp,1))
			c = model.comp[idx,2:11]
			#println(c)
			s = model.seismic[idx,2:end]
			p = pdf(model.er, s .- [rho[samp],vp[samp],vpvs[samp]])
			p = p*pdf(model.currentPrior, c)/pdf(model.basePrior, c)
			maxp = pdf(model.er, model.er.μ) # make alg more efficient by limiting rand vals from 0 to max p
			# if isnan(p)
			# 	println("prob is nan, comp $c, seismic props $(s)")
			# end
			#println("prob = $(p*100000)/100000, maxp = $maxp")
			if rand()*maxp < p
				toreturn[samp,:] .= model.comp[idx,:]
				#println("Selected sample $(model.comp[idx])")
				break
			end
			if i == 100000
				println("Did not find sample for seismic props $([rho[samp],vp[samp],vpvs[samp]])")
			end
		end
	end

	return toreturn, zeros(size(toreturn))
end

# Just returns one elemnt of interest
function resultSize(model::RejectionPriorModel)
	return length(PERPLEX_ELEMENTS) ##
end

"""
Output relevant model info.
"""
function modelSummary(model::RejectionPriorModel)
	return "Means of priors\n"*join(model.means, "\n")
end
