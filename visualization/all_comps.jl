"""
For a final result file, plot composition over crust base depth, binned by age. 
"""

using ArgParse
using DelimitedFiles 
using Plots; gr();
using Statistics

include("../src/config.jl")
include("../bin.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/latlong_weighted"
    "--compare"
    	help = "Look for rf/crust1 comparison results instead of range model results"
    	arg_type = Bool 
    	default = false 
end 
parsed_args = parse_args(ARGS, s)

if parsed_args["compare"] == false 
	files = ["data/"*parsed_args["data_prefix"]*"/output/results-"*layer*"-range-earthchem.csv" for layer in LAYER_NAMES]
else 
	files = vcat(["data/"*parsed_args["data_prefix"]*"/output/results-rf_vs_crust1-crust1-$(layer).csv" for layer in LAYER_NAMES],
		["data/"*parsed_args["data_prefix"]*"/output/results-rf_vs_crust1-rf-$(layer).csv" for layer in LAYER_NAMES])
end
outputPath = "data/"*parsed_args["data_prefix"]*"/output"

for (l, file) in enumerate(files)	
	if l%3 == 0 
		layer = LAYER_NAMES[3]
	else 
		layer = LAYER_NAMES[l%3]
	end 

	(results, header) = readdlm(file, ',', header=true)
	header = header[:] # make 1d 

	si_index = findfirst(isequal("SiO2"), header) # result 
	co_index = findfirst(isequal("CO2"), header)

	good = .~ isnan.(sum(results[:,si_index:co_index], dims=2))
	good = good[:]
	means = median(results[good,si_index:co_index], dims=1)
	println("For layer $layer, $(size(results,1) - sum(good)) test points have at least one NaN result value")

	println("$layer")
	println("$(header[si_index:co_index])")
	println(join(round.(means, sigdigits=3), " & "))


	# 90th and 10th quantiles 
	res = zeros((2,0))
	for c in eachcol(results[good,si_index:co_index])
		res = hcat(res, quantile(c, [.1,.9]))
	end 

	println("10th and 90th percentile cutoffs:")
	println(join(round.(res[1,:], sigdigits=3), ", "))
	println(join(round.(res[2,:], sigdigits=3), ", "))

end 