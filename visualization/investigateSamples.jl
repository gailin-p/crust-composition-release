# For a dataset from commutativityOfAveraging

using ArgParse
using JLD
using MAT
using HDF5
using Statistics
using StatGeochem
using MultivariateStats
using Interpolations
using StatsBase
using Plots; gr();

include("../bin.jl") # TODO locate and use statgeochem version of this
include("../src/utilities.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data"
        help = "JLD data file saved by commutativityOfAveraging"
        arg_type = String
        required = true
    "--make_figs"
        help = "Figures?"
        arg_type = Bool
        default = false
    "--prefix"
        help = "Prefix for output files, if any (if make_figs or mcign/pca provided)"
        arg_type = String
        default = "../output/average/"
end
parsed_args = parse_args(ARGS, s)

# order of data: layer, property, run
layers = ["Upper","Middle","Lower"]
props = ["rho","vp","vpvs"] # properties are in this order.

# load data
data = load(parsed_args["data"])
props_of_ave = data["props_of_ave"]
ave_props = data["ave_properties"]
indices = data["indices"]
r = size(indices)[2] # num samples
n = size(indices)[1] # samples per average

# Compare averages
for (l, layer) in enumerate(layers)
    for (p, prop) in enumerate(props)
        ave1 = inverseMean(props_of_ave[l,p,:])
        sigma1 = nanstd(props_of_ave[l,p,:])/sqrt(r)
        ave2 = inverseMean(ave_props[l,p,:])
        sigma2 = nanstd(ave_props[l,p,:])/sqrt(r)
        Z = (ave1-ave2)/sqrt(sigma1+sigma2)
        println("Z = $Z for $layer, $prop")
    end
end

# Compare averages
for (l, layer) in enumerate(layers)
    for (p, prop) in enumerate(props)
        dif = props_of_ave[l,p,:] .- ave_props[l,p,:]
        percent = 100 .* abs.(dif ./ props_of_ave[l,p,:])
        within5 = sum(percent .< 5)
        println("For $layer, $prop : 
            Average difference = $(nanmean(dif)), with std $(nanstd(dif))
            Average % dif $(nanmean(percent)) and max % dif $(maximum(percent))
            $(100*within5/length(percent)) % are less than 5% difference.")
    end
end

if parsed_args["make_figs"]
    diffs = props_of_ave - ave_props # layer, property, run (3,3,r)

    for (l, layer) in enumerate(layers)
        for (p, prop) in enumerate(props)
            left = min(minimum(props_of_ave[l,p,:]), minimum(ave_props[l,p,:]))
            right = max(maximum(props_of_ave[l,p,:]), maximum(ave_props[l,p,:]))

            # Stacked histograms: ave of props, props of ave, diffs
            h1 = histogram(props_of_ave[l,p,:], legend=false, xlabel="Property of average", bins=60, xlims=(left,right)) # xlims
            h2 = histogram(ave_props[l,p,:], legend=false, xlabel="Average properties", bins=60, xlims=(left,right))
            h3 = histogram(diffs[l,p,:], legend=false, xlabel="$(prop) difference", bins=60)

            h = plot(h1, h2, h3, layout=(3,1), size=(800,600))

            savefig(h, parsed_args["prefix"]*"compare_hist_$(prop)_$(layer)_n$(n)_r$(r).pdf")
        end
    end
end
