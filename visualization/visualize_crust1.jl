"""
Make analagous pointclouds to visualize_range_model, so we can see where earthchem data falls compared to our models. 

Then run model on data and output same plots but colored

TODO color points by geotherm 
"""

using Plots; gr();
using ArgParse 

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="test"
    "--model", "-m"
        help = "Type of model to use. Allowed: inversion, range"
        arg_type = String
        range_tester = x -> (x in ["inversion","range"])
        default = "inversion"
end
parsed_args = parse_args(ARGS, s) # TODO currently unused except for deciding where to place output 

include("../src/config.jl")
include("../src/crustDistribution.jl")
include("../src/inversionModel.jl")

outputPath = "data/"*parsed_args["data_prefix"]*"/output/not_model_specific"
mkpath(outputPath)

for (i, layer) in enumerate(LAYER_NAMES)
	data = crustDistribution.getAllSeismic(i+5, resample=false)
	data = unique(hcat(data[1:end-1]...), dims=1) # only 20 unique combos in unresampled data 

	# 3d scatter plot 
	cameras = ((30,30), (30,80), (80,10))
	plots = []
	for (i, camera) in enumerate(cameras)
		push!(plots, plot(data[:,1], data[:,2], data[:,3], camera=camera,
			seriestype=:scatter, markersize=5, markerstrokewidth=0, legend=false,
			colorbar=false)) # labels (xlabel="Rho (kg/m^3)", ylabel="Vp (m/s)", zlabel="Vp/Vs") look awful
	end 
	p = plot(plots..., size=(600,1500), layout=(3,1))
	savefig(p, "$outputPath/crust1_3dplot_$layer.png")
end

# Now run model and color according to model results 
modelOutputPath = "data/"*parsed_args["data_prefix"]*"/output/$(parsed_args["model"])_model"
mkpath(modelOutputPath)

if parsed_args["model"] == "range"
	models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
elseif parsed_args["model"] == "inversion"
	models = makeModels(parsed_args["data_prefix"], modelType=InversionModel)
end 

for (i, layer) in enumerate(LAYER_NAMES)
	data = crustDistribution.getAllSeismic(i+5, resample=true, n=100)

	results, _ = estimateComposition(models, layer, data...)

	# 3d scatter plot 
	cameras = ((30,30), (30,80), (80,10))
	plots = []
	for (i, camera) in enumerate(cameras)
		colorbar=false
		if i == 3 # colorbar on last plot only
			colorbar=:bottom 
		end
		push!(plots, plot(data[1:end-1]..., camera=camera,
			marker_z=results, seriestype=:scatter, markersize=1, markerstrokewidth=0, legend=false,
			colorbar=colorbar, colorbar_title="SiO2 (wt %)")) # labels (xlabel="Rho (kg/m^3)", ylabel="Vp (m/s)", zlabel="Vp/Vs") look awful
	end 
	p = plot(plots..., size=(600,1500), layout=(3,1))
	savefig(p, "$modelOutputPath/testPoints_results_$layer.png")
end



