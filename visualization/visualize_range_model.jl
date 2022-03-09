"""
For binned models, visualizes 1/nbins samples from each bin.
"""

using ArgParse
using Plots; gr();
using Statistics
using DelimitedFiles

include("../src/config.jl")
include("../src/inversionModel.jl")
include("../src/invertData.jl")
include("../src/crustDistribution.jl")
include("../src/rejectionModel.jl")
include("../src/cracks.jl")

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
    	default=false
    "--show_inverted"
    	help = "Plot actual inverted data points on top of model plot?"
    	arg_type=Bool
    	default=false
end
parsed_args = parse_args(ARGS, s)


# crackFile = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/crack_profile.csv"
# if !isfile(crackFile)
# 	crackFile = ""
# end

option_file = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/inversion_options.csv"
p = Dict()
params = readdlm(option_file, ',', header = false)
for row in 1:size(params,1)
    p[params[row,1]] = params[row,2]
end
println(p)
crack_porosity = p["fraction_crack"] * p["crack"]
pore_porosity = (1-p["fraction_crack"]) * p["crack"]
profile = CrackProfile(sphere_properties,
	(p["wet_crack"] ? water_rho : 0), #
	(p["wet_crack"] ? water_K : 0), # 0 for dry cracks
	pore_porosity,
	crack_porosity,
	(p["wet_crack"] ? "wet" : "dry"),
	"sphere",
	p["alteration_fraction"] # fraction of clay alteration minerals)
)

models = makeModels(parsed_args["data_prefix"], modelType=RejectionModel, crackProfile=profile)

outputPath = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/output"
mkpath(outputPath)

layer_plots = []
for (l, layer) in enumerate(LAYER_NAMES)
	# 3d scatter plot
	#model = models.models[layer][(max(1,floor(Int,models.nbins/2)))] # use an intermediate geotherm
	cameras = ((30,30), (30,80), (80,10))
	plots = []

	## bug in gr makes gif(.) function hang in some julia contexts, so doing this in a notebook instead
	# Make gif of rotating camera
	# model = models.models[layer][3] # use intermediate geotherm for viz purposes
	# subsample = sample(1:size(model.seismic,1), 1000) # indices
	# anim = @animate for i = 10:70
	# 	println("anim $i")
	# 	plot(model.seismic[subsample,2], # rho
	# 		model.seismic[subsample,3], model.seismic[subsample,4], # vp, vp/vs
	# 		camera=(i,60), marker_z=model.comp[subsample,2], zlims=(1.25,2.0),
	# 		xlims=(2250,4000), ylims=(4.75,8.25),
	# 		seriestype=:scatter, markersize=2, markerstrokewidth=0)
	# end
	# gif(anim, "$outputPath/$(layer).gif", fps = 10)

	# Make still imgs at various angles
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
			seriestype=:scatter, markersize=5, markeralpha=1, markerstrokewidth=0, legend=false, markercolor=:green2)
		end
		if parsed_args["show_inverted"] # Plot inverted points on same axes
			resdat, resh = readdlm("data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/results-$layer.csv", ',', header=true)
			plot!(q, resdat[:,1], resdat[:,2], resdat[:,3],
			seriestype=:scatter, markersize=2, markeralpha=1, markerstrokewidth=0, legend=false, markercolor=:cyan)
		end
		for model in [models.models[layer][2]]
			total = size(model.seismic,1)
			plot!(q,
				#model.seismic[1:models.nbins:total,2], # rho
				#model.seismic[1:models.nbins:total,3], #model.seismic[1:models.nbins:total,4], # vp, vp/vs
				model.seismic[:,2], # rho
				model.seismic[:,3],
				model.seismic[:,4], # vp, vp/vs
				camera=camera,
				#marker_z=model.comp[1:models.nbins:total,2],
				marker_z=model.comp[:,2], 
				zlims=(1.25,2.0),
				xlims=(2250,4000), ylims=(4.75,8.25),
				seriestype=:scatter, markersize=2, markerstrokewidth=0) # labels (xlabel="Rho (kg/m^3)", ylabel="Vp (m/s)", zlabel="Vp/Vs") look awful
		end
		if i==1
			push!(layer_plots, q)
		end
		push!(plots, q)
	end
	p = plot(plots..., size=(600,1500), layout=(3,1))
	savefig(p, "$outputPath/3dplot_bin2_$layer.png")
end

r = plot(layer_plots..., size=(1500, 600), layout=(1,3))
savefig(r, "$outputPath/3dplot_all_layers_bin2.png")
