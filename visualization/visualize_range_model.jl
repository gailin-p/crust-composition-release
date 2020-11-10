"""
For binned models, visualizes 1/nbins samples from each bin. 
"""

using ArgParse
using Plots; gr();

include("../src/inversionModel.jl")
include("../src/invertData.jl")
include("../src/config.jl")
include("../src/crustDistribution.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for perplex output"
        arg_type= String
        default="test"
    "--name"
        help = "Name of model run"
        arg_type= String
        default="default"
    "--show_test"
    	help = "Plot Crust1.0 test points on top of model plot?"
    	arg_type=Bool
    	default=true
end
parsed_args = parse_args(ARGS, s)


crackFile = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/crack_profile.csv"
if !isfile(crackFile)
	crackFile = ""
end
models = makeModels(parsed_args["data_prefix"], modelType=RangeModel, crackFile=crackFile)

# histograms of vp, vs, rho for each model 
plots = [] 
for layer in LAYER_NAMES 
	for model in models.models[layer]
		push!(plots, histogram(model.seismic[:,2], legend=false, xlims=(2000,4000)));
	end 
end 

outputPath = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/output/range_model"
mkpath(outputPath)

# for (i, p) in enumerate(plots)
# 	savefig(p, "$outputPath/rho-hist-$i.pdf")
# end

p = plot(plots..., size=(1000,2000), layout=(models.nbins, length(LAYER_NAMES))); 

savefig(p, "$outputPath/rho-hist.pdf")

layer_plots = []
for (l, layer) in enumerate(LAYER_NAMES)
	# 3d scatter plot 
	#model = models.models[layer][(max(1,floor(Int,models.nbins/2)))] # use an intermediate geotherm 
	cameras = ((30,30), (30,80), (80,10))
	plots = []
	for (i, camera) in enumerate(cameras)
		colorbar=false
		# if (i == 3) | (l==3) # colorbar on last plot only (for all layers, on bottom layer.)
		# 	colorbar=:bottom 
		# end
		q = plot(legend=false, colorbar=colorbar)
		# Add 1/nbins of points from every model to get overall view 
		if parsed_args["show_test"] # Plot Crust1.0 test points on same axes 
			tests, _ = getAllSeismic(l+5, resample=false)
			tests = unique(hcat(tests[1:end-1]...), dims=1) # only 20 unique combos in unresampled data 
			plot!(q, tests[:,1], tests[:,2], tests[:,3], 
			seriestype=:scatter, markersize=8, markeralpha=1, markerstrokewidth=0, legend=false, markercolor=:green2)
		end
		for model in models.models[layer]
			total = size(model.seismic,1)
			plot!(q, model.seismic[1:models.nbins:total,2], # rho
				model.seismic[1:models.nbins:total,3], model.seismic[1:models.nbins:total,4], # vp, vp/vs 
				camera=camera, marker_z=model.comp[1:models.nbins:total,2], zlims=(1.5,3.5),
				seriestype=:scatter, markersize=1, markerstrokewidth=0) # labels (xlabel="Rho (kg/m^3)", ylabel="Vp (m/s)", zlabel="Vp/Vs") look awful
		end
		if i==1 
			push!(layer_plots, q)
		end
		push!(plots, q)
	end 
	p = plot(plots..., size=(600,1500), layout=(3,1))
	savefig(p, "$outputPath/3dplot_$layer.png")
end

r = plot(layer_plots..., size=(1500, 600), layout=(1,3))
savefig(r, "$outputPath/3dplot_all_layers.png")
