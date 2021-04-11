"""
Central script for algorithm 
"""

using DelimitedFiles
using Distributions
using GaussianMixtures
using MultivariateStats
using ProgressMeter
using StatGeochem

include("../src/config.jl")
include("../src/crustDistribution.jl")
include("../src/seismic.jl")
include("../src/parsePerplex.jl")

struct PerplexConfig
	solutions::String
	dataset::String
	perplex::String
	scratch::String
	fluid_endmembers::String 
	npoints::Integer
	st::SeismicTransform 
end

struct MCMCRunner
	std::Array{Float64,1}
	prior::MixtureModel
	config::PerplexConfig
	er::MvNormal
end

struct DataSample
	tc1::Float64 # depth to 550 
	depth::Float64 # depth of sample / middle of layer 
	vp::Float64 # target seismic properties: vp, vpvs, rho 
	vpvs::Float64
	rho::Float64
end

function run_mcmc(n::Integer, runner::MCMCRunner, data::DataSample) 
	m = length(COMPOSITION_ELEMENTS)
	samples = fill(NaN, (n, m)) 
	accepted = fill(false, n) # for calculating acceptance ratio and seeing burn-in
	logprob = fill(NaN, n)

	# Choose first sample and get its P(x)
	# First samples is the mean of lat-long resampled earthchem 
	x = normalize([58.67910944210001, 0.9524187580000001, 14.5812233658, 7.2417907846, 4.8595292555, 5.957280192800002, 3.183153271, 2.1429591833, 1.8062233445, 0.5963123932000001])
	Px = probability(runner,data, x)
	if isnan(Px)
		Px = -Inf
	end 
	println("Comparing all to $Px")

	# loop 
	@showprogress 1 "Iterating..." for i in 1:n
		y = proposal(runner, x) # proposal 
		Py = probability(runner, data, y)
		
		# Accept or reject 
		if (Py > Px) | (log(rand()) <= (Py - Px)) # accept 
			x .= y 
			Px = Py 
			accepted[i] = true 
		end

		# save sample 
		logprob[i] = Px
		samples[i,:] .= x
	end

	return samples, accepted, logprob
end

function PerplexConfig(perplex, scratch)
	solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)"*
		"\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar"*
		"\nDo(HP)\n"
	npoints = 20
	fluid_endmembers = "abL\nanL\ndiL\nenL\nfaL\nfliq\nfoL\nkspL\nmliq\nqL\nsiL\nq8L\nfa8L\nfo8L\nsil8L\nh2oL\nh2o8L\n"
	
	dataset = "hpha11ver.dat"

	st = SeismicTransform() 

	return PerplexConfig(solutions, dataset, perplex, scratch, fluid_endmembers, npoints, st)
end 


function PerplexConfig()
	return PerplexConfig("/Users/gailin/resources/perplex-stable/", "/Users/gailin/dartmouth/crustal_structure/perplexed_pasta/")
end

"""
Use some local default locations of resampled earthchem and perplex error files 
"""
function MCMCRunner(perplex_config=PerplexConfig())
	# Prior 
	ign, h = readdlm("../data/remote/base_nobin/bsr_ignmajors_1.csv", ',', header=true)
	ign = ign[:,2:11] # no index, no context 
	n = 1000 # does not converge for all 100000 samples
	gm = GMM(11, Matrix(ign[rand(1:size(ign)[1],n),:]))
	mw = MixtureModel(gm);

	# error distribution 
	dat, header = readdlm("../data/adjustedPerplexRes.csv", ',', header=true);
	header = header[:]
	# Interested in non-nan differences between actual and Perplex data 
	ers = fill(NaN, (size(dat)[1],3));
	ers[:,1] .= dat[:,findfirst(isequal("perplex vp"), header)] .- dat[:,findfirst(isequal("dabie vp"), header)]
	ers[:,2] .= dat[:,findfirst(isequal("perplex vp/vs"), header)] .- dat[:,findfirst(isequal("dabie vp/vs"), header)]
	ers[:,3] .= dat[:,findfirst(isequal("perplex rho"), header)] .- (1000 .* dat[:,findfirst(isequal("dabie rho"), header)]);
	ers = ers[.!(isnan.(sum(ers, dims=2)))[:],:]; 
	N_er = fit(MvNormal, ers')

	# updates as fraction of avg composition 
	update = [58.67910944210001, 0.9524187580000001, 14.5812233658, 7.2417907846, 4.8595292555, 5.957280192800002, 3.183153271, 2.1429591833, 1.8062233445, 0.5963123932000001] ./ 10

	return MCMCRunner(update, mw, perplex_config, N_er)
end

function normalize(elts::Array{Float64,1})
	elts[elts .< 0] .= 0.0
	return 100 .* elts ./ sum(elts)
end

function run_perplex(config::PerplexConfig, data::DataSample, x::Array{Float64,1})
# run perplex on x 
	depth = data.depth # middle of layer 
	dtdz = 550.0/data.tc1 # sampling geotherm 
	formation_depth, formation_dtdz = crustDistribution.getFormationParams(depth, dtdz)

	geotherm = formation_dtdz/dpdz # dt/dp 
	P_range = [formation_depth*(9/10)*dpdz, formation_depth*(11/10)*dpdz] # we only need a small range around formation t, p

	# Run perplex
	perplex_configure_geotherm(config.perplex, config.scratch, x, PERPLEX_COMPOSITION_ELTS,
        P_range, 273.15, geotherm, dataset=config.dataset, solution_phases=config.solutions,
        excludes=config.fluid_endmembers, index=1, npoints=config.npoints)
    point = perplex_query_point(config.perplex, config.scratch, formation_depth*dpdz, index=1)
    #println("input conditions depth = $depth, temp = $(depth*dtdz) C")
    #println("conditions depth = $formation_depth, temp = $(formation_depth*formation_dtdz) C")
    #println(point)

    try 
		properties = get_system_props(point)
		endmembers = parse_perplex_point(point)

		P = dpdz*depth + 280 # add surface pressure (bar)
		T = dtdz*depth + 273.15 # add surface temp (K)
		rho, vp, vs = get_seismic(T, P, properties, endmembers, config.st)
		return [vp, vp/vs, rho]  
	catch e
		println("\r\n\r\nCannot process sample due to \r\n $e")
		return [NaN, NaN, NaN]
	end 
end 

# for acceptance ratio 
function probability(runner::MCMCRunner, data::DataSample, x::Array{Float64,1})
	prior = logpdf(runner.prior, x)

	props = run_perplex(runner.config, data, x)
	update = logpdf(runner.er, props .- [data.vp, data.vpvs, data.rho])

	return prior + update 
end


function proposal(runner::MCMCRunner, x::Array{Float64,1})
	return normalize((randn(length(x)) .* runner.std) .+ x)
end



