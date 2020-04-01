"""
For a final result file, plot composition over crust base depth, binned by age. 
"""

using ArgParse
using DelimitedFiles 
using Plots; gr();

include("../src/config.jl")
include("../bin.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/latlong_weighted"
end 
parsed_args = parse_args(ARGS, s)

files = ["data/"*parsed_args["data_prefix"]*"/output/results-"*layer*"-range-earthchem.csv" for layer in LAYER_NAMES]
outputPath = "data/"*parsed_args["data_prefix"]*"/output"

for (l, file) in enumerate(files)	
	(results, header) = readdlm(file, ',', header=true)
	header = header[:] # make 1d 

	si_index = findfirst(isequal("SiO2"), header) # result 
	age_index = findfirst(isequal("sample_age"), header) # of input data (crust1.0 data point)
	depth_index = findfirst(isequal("sample_depth"), header) # of crust1.0 data point 
		
	minBase = minimum(results[:,depth_index])
	maxBase = maximum(results[:,depth_index])

	ageBins = [-50,50,2500,4500]

	if l>1
		p = plot(xlabel="Base of Crust1.0 (km)", ylabel="Estimated silica (wt%)", legend=:bottomright)
	else 
		p = plot(xlabel="Base of Crust1.0 (km)", ylabel="Estimated silica (wt%)")
	end

	for (i, base) in enumerate(ageBins[1:end-1])
		top = ageBins[i+1]
		test = (results[:,age_index] .> base) .& (results[:,age_index] .<= top)
		centers, aves, yerrors = bin(results[test, depth_index], results[test,si_index], minBase, maxBase, 10) # x,y,min,max,nbins
		plot!(p, centers, aves, yerror=yerrors, label="$(LAYER_NAMES[l]) $(base)Ma - $(top)Ma", markerstrokecolor=:auto)
	end 

	savefig(p, outputPath*"/composition-crustbase-$(LAYER_NAMES[l])-ageBinned.pdf")

end 















