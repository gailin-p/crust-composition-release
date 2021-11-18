"""
Create main figure (composition of upper, middle, and lower crust over time)

"""

using ArgParse
using DelimitedFiles
using Plots; gr();
using Statistics
using StatsBase

include("../src/bin.jl")
include("../src/config.jl")
include("../src/crustDistribution.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_path", "-p"
        help = "Path to model output"
        arg_type= String
        required=true
    "--include_prior"
        help = "Include histogram of prior means?"
        arg_type = Bool
        default= true
end

parsed_args = parse_args(ARGS, s)


build_args = readdlm("$(parsed_args["data_path"])/inversion_options.csv", ',', header=false)
N = build_args[findfirst(isequal("num_invert"), build_args[:,1]),2]
M = build_args[findfirst(isequal("num_runs"), build_args[:,1]),2]

outputPath = "$(parsed_args["data_path"])/output"
mkpath(outputPath) # make if does not exist

# Sample M earths of size N from prior compositions: this would be the estimated earth comp from the prior with no info from seismic properties
if parsed_args["include_prior"]
    ign, h = readdlm("$(parsed_args["data_path"][1:findlast("/",parsed_args["data_path"])[1]])bsr_ignmajors_1.csv", ',', header=true)
    si_idx = findfirst(isequal("SiO2"), h[:])
    prior_earths = zeros(M)
    for i in 1:M
        samples = sample(1:size(ign,1), N)
        prior_earths[i] = mean(ign[samples,si_idx])
    end
end

## Compare to R&G
ures, h = readdlm("$(parsed_args["data_path"])/results-upper.csv", ',', header=true)
mres, h = readdlm("$(parsed_args["data_path"])/results-middle.csv", ',', header=true)
lres, h = readdlm("$(parsed_args["data_path"])/results-lower.csv", ',', header=true)
si_i = findfirst(isequal("SiO2"), h[:]);

resu = zeros(M,3)

for (l, res) in enumerate([ures, mres, lres])
    for i in 1:M
        istart = (i-1)*N + 1
        iend = i*N
        #println(nanmean(res[istart:iend, si_i]))
        resu[i,l] = nanmean(
            convert(Array{Float64,1}, res[istart:iend, si_i]))
    end
end

colors = [:blue, :orange, :green]
stephist(resu[:,1], normalize=:pdf, label="Upper", c=colors[1], nbins=15)
stephist!(resu[:,2], normalize=:pdf, label="Middle", c=colors[2], nbins=15)
stephist!(resu[:,3], normalize=:pdf, label="Lower", c=colors[3], nbins=15)
if parsed_args["include_prior"]
    stephist!(prior_earths, normalize=:pdf, label="Prior", nbins=15)
end
scatter!([66.6], [1/100], c=:black, shape=:diamond, label="Rudnick & Gao (2014)", legend=:outerright)
scatter!([66.6], [1/100], c=colors[1], label=false, shape=:diamond)
scatter!([63.5], [1/100], c=colors[2], label=false, shape=:diamond)
scatter!([53.4], [1/100], c=colors[3], label=false, shape=:diamond)
plot!(yticks=false, framestyle=:box, xlabel="SiO2", size=(600, 200),
    title="$(parsed_args["data_path"])", xlims=(47.5, 67.5))
savefig("$(parsed_args["data_path"])/output/compare_rg.pdf")
