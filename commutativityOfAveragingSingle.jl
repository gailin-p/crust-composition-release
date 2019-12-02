### Single threaded version of commutativityOfAveragingMPI

using ArgParse
using Plots; gr();
using JLD

include("src/commutativityOfAverage.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/scratch/gailin/" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        required = true
    "--ign", "-i"
        help = "Path to ignmajors file"
        arg_type = String
        default = "ignmajors.csv"
    "--n_samples", "-n"
        help = "Number of samples to average"
        arg_type = Int
        default = 2
    "--n_runs", "-r"
        help = "How many sample pairs to run"
        arg_type = Int
        default = 30
end
parsed_args = parse_args(ARGS, s)
r_requested = parsed_args["n_runs"]
n = parsed_args["n_samples"]
r = parsed_args["n_runs"]

ave_properties = Array{Float64,3}(undef,(3,3,r))
props_of_ave = Array{Float64,3}(undef,(3,3,r))
indices = Array{Int,2}(undef,(n,r)) # store n indices for each of r runs, for reconstruction!

# run samples
commutativityOfAverage.runSamples!(props_of_ave, ave_properties, indices,
            r, n,
            parsed_args["ign"], parsed_args["perplex"], parsed_args["scratch"])

# Save data and plots
# Plot diffs
diffs = props_of_ave-ave_properties
layers = ["upper","middle","lower"]
props = ["rho","vp","vpvs"]
for l in 1:3
    for p in 1:3
       h = histogram(diffs[l,p,:], legend=false, xlabel="$(props[p]) difference")
       savefig(h, "output/average/diff_hist_$(props[p])_$(layers[l])_n$(n)_r$(r).pdf")
    end
end

save("data/commutativity-n$(n)-r$(r).jld","ave_properties", ave_properties, "props_of_ave", props_of_ave, "indices", indices)
