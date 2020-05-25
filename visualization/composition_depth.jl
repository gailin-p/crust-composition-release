"""
Is the change in composition of the lower crust over age a result of changing crust thickness over age? 

Get (binned) si/base of crust relationship for lower crust using modern ()

Use si/base relationship to predict crust si from crust depth (from Crust1.0)

Plot si predicted from si/base relationship vs si predicted by our model. 
"""

using ArgParse
using DelimitedFiles 
using Plots; gr();
using Interpolations

include("../src/config.jl")
include("../src/bin.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/latlong_weighted"
    "--model_run"
    	help="Which model run suffix? eg -range-earthchem for range model, earthchem age model"
    	arg_type= String 
    	default="-range-earthchem"
end 
parsed_args = parse_args(ARGS, s)

files = ["data/"*parsed_args["data_prefix"]*"/output/results-"*layer*parsed_args["model_run"]*".csv" for layer in LAYER_NAMES]
outputPath = "data/"*parsed_args["data_prefix"]*"/output"

l = 3 # lower 
file = files[l]
(results, header) = readdlm(file, ',', header=true)
header = header[:] # make 1d 

si_index = findfirst(isequal("SiO2"), header) # result 
age_index = findfirst(isequal("sample_age"), header) # of input data (crust1.0 data point)
depth_index = findfirst(isequal("sample_depth"), header) # of crust1.0 data point 
	
# si/base relationship 
minBase = minimum(results[:,depth_index])
maxBase = maximum(results[:,depth_index])
ageBin = [-50,50] # consider these ages for si/base binning 
test = (results[:,age_index] .> ageBin[1]) .& (results[:,age_index] .<= ageBin[2])
c, m, yerrors = bin(results[test, depth_index], results[test,si_index], minBase, maxBase, 20) # x,y,min,max,nbins

# Interpolation of si/base relationship
itp = interpolate(m, BSpline(Linear())) # Linear interpolation of bin means
itp = scale(itp, c) # Scale by bin centers
itp = extrapolate(itp, Flat()) # outside bin means, just give the upper or lower bin mean. 

# For all samples older than 50ma, get base and si. 
age_max = 4000
test_old = (results[:,age_index] .> ageBin[2]) .& (results[:,age_index] .<= age_max)
depths = results[test_old, depth_index]
sis = results[test_old, si_index]
depth_predicted_si = itp(depths)

scatter(depth_predicted_si, sis, marker_z=depths, label="Samples", markersize=2, markerstrokewidth=0,
	ylabel="Si predicted by model", xlabel="Si predicted from crust base depth", colorbar_title="Sample depth")

savefig(outputPath*"/depth_si_vs_model_si.pdf")
