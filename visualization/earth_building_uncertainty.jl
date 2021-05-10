"""
Takes N random samples the size of the continents (assuming each sample covers 1x1 deg, ie Earthchem)
Finds the percentiles of the compositions of those samples 

"""

using ArgParse
using DelimitedFiles 
using Plots; gr();
using Statistics
using ProgressMeter
using StatsBase

include("../src/bin.jl")
include("../src/config.jl")
include("../src/crustDistribution.jl")
include("../src/invertData.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_path", "-p"
        help = "Path to model output"
        arg_type= String
        required=true
    "--n_globes", "-n"
        help = "How many globes to make?"
        arg_type = Int
        default = 10000
end 

parsed_args = parse_args(ARGS, s)

println("Choosing $(parsed_args["n_globes"]) globes x 3 layers of size $(nOriginal())")


# By age! 
function make_globes(file, N)
	globe_si = Array{Float64, 1}(undef, N) 

	(result, header) = readdlm(file, ',', header=true)
	header = header[:] # make 1d 

	si_index = findfirst(isequal("SiO2"), header) # result 
	#good = .~isnan.(result[:,si_index])

	globe_n = nOriginal() ## Number of 1x1 deg squares in continents 
	sample_p = globe_n/size(result,1)
	for g in 1:N 
		indices = randsubseq(1:size(result,1), sample_p)
		#println("chose $(length(indices)) indices")
		globe_si[g] = nanmean(Float64.(result[indices, si_index]))
	end

	return globe_si
end 

p = []

files = ["$(parsed_args["data_path"])/results-$layer.csv" for layer in LAYER_NAMES]
for (l, file) in enumerate(files)
	globes = make_globes(file, parsed_args["n_globes"])
	println(mean(globes))
	println("95 percentiles $(percentile(globes, 5)) - $(percentile(globes, 95))")
	push!(p, histogram(globes, title=LAYER_NAMES[l], legend=false, normalize=:pdf, xlims=(53,67)))
end

p1 = plot(p..., layout=(3,1), size = (400,600))
outputPath = "$(parsed_args["data_path"])/output"
mkpath(outputPath) # make if does not exist
savefig(p1, outputPath*"/world-building.pdf");


# Upper, middle, and lower from Rudnick and Fountain 
# plot!(p, [200], [66.0], seriestype=:scatter, marker=8,markershape=:star5, markercolor=:blue, label="")
# plot!(p, [200], [60.6], seriestype=:scatter, marker=8,markershape=:star5, markercolor=:orange, label="")
# plot!(p, [200], [52.3], seriestype=:scatter, marker=8,markershape=:star5, markercolor=:green, label="")


# Plot exposed. Always use base so have same comparison... 
# TODO this will need to change when re-do base w new elts
# tmin = 0
# tmax = 4000
# nbins = 4
# samplePath = "data/remote/base/bsr_ignmajors_1.csv"
# ign, h = readdlm(samplePath, ',', header=true)
# si_i = findfirst(isequal("SiO2"),PERPLEX_ELEMENTS)
# mean_exposed_si = nanmean(ign[:,si_i])





