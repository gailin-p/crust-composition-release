"""
For a final result file, print composition means
"""

using ArgParse
using DelimitedFiles
using Plots; gr();
using Statistics
using StatsBase

include("../src/config.jl")
include("../src/bin.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_path", "-p"
        help = "Folder for data files"
        arg_type= String
        default="remote/base/default"
    "--compare_rg"
    	help = "Compare to Rudnick and Gao"
    	arg_type = Bool
    	default = false
end
parsed_args = parse_args(ARGS, s)

build_args = readdlm("$(parsed_args["data_path"])/inversion_options.csv", ',', header=false)
N = build_args[findfirst(isequal("num_invert"), build_args[:,1]),2]
M = build_args[findfirst(isequal("num_runs"), build_args[:,1]),2]

#         SiO2 TiO2 Al2O3 FeO  MgO  CaO  Na2O K2O
rg_dat = [66.6 0.64 15.4  5.04 2.48 3.59 3.27 2.80
		  63.5 0.69 15.0  6.02 3.59 5.25 3.39 2.30
		  53.5 0.82 16.9  8.57 7.24 9.59 2.65 0.61]


files = ["$(parsed_args["data_path"])/results-$layer.csv" for layer in LAYER_NAMES]

if parsed_args["compare_rg"] == true
	this_dat = zeros((3,8))
end

for (l, file) in enumerate(files)
	if l%3 == 0
		layer = LAYER_NAMES[3]
	else
		layer = LAYER_NAMES[l%3]
	end

	if l%3 == 1
		#println(file)
	end

	(results, header) = readdlm(file, ',', header=true)
	header = header[:] # make 1d

	si_index = findfirst(isequal("SiO2"), header) # result
	co_index = findfirst(isequal("CO2"), header)

	good = .~ isnan.(sum(results[:,si_index:co_index], dims=2))
	good = good[:]
	means = mean(results[good,si_index:co_index], dims=1)

	runs_means = zeros((M,10))
	for m in 1:M
		start = (m-1)*N + 1
		eend = m*N
		result_subset = results[start:eend,si_index:co_index]
		good_subset = .~ isnan.(sum(result_subset, dims=2))
		good_subset = good_subset[:]
		runs_means[m,:] .= mean(result_subset[good_subset,:], dims=1)[:]
	end

	percentile_low = zeros(10)
	percentile_high = zeros(10)
	for i in 1:10
		percentile_low[i] = percentile(runs_means[:,i], 5)
		percentile_high[i] = percentile(runs_means[:,i], 95)
	end
	#println("For layer $layer, $(size(results,1) - sum(good)) test points have at least one NaN result value")

	println("$layer")
	#println("$(header[si_index:co_index])")
	println(join(round.(means, sigdigits=3), " & "))
	println(join(round.(percentile_low, sigdigits=3), " & "))
	println(join(round.(percentile_high, sigdigits=3), " & "))

	if parsed_args["compare_rg"]
		this_dat[l,:] = means[1:end-2]
	end

end

if parsed_args["compare_rg"]
	println("Difference as % of this result")

	dif = 100 .* (this_dat .- rg_dat)./this_dat
	for row in 1:size(dif,1)
		println(join(round.(dif[row,:], sigdigits=3), " & "))
	end

	println("Difference as wt %")

	dif = this_dat .- rg_dat
	for row in 1:size(dif,1)
		println(join(round.(dif[row,:], sigdigits=3), " & "))
	end
end
