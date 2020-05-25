"""
Compare our results to results of other workers 

TODO: add standard error of mean bars. 
"""

using Plots; gr();
using Statistics 
include("../src/crustDistribution.jl")

mutable struct Result
	name::String
	composition::Array{Float64, 1} # Lower, middle, upper
	layers::Array{Float64, 1} # Depth or layer index?
	#depth::Array{Float64, 1} # Depths for estimates 
	notes::String
end 

mutable struct ResultWithError
	name::String
	composition::Array{Float64, 1} # Lower, middle, upper
	sem::Array{Float64, 1} # standard error of mean; lower, middle, upper
	layers::Array{Float64, 1} # Depth or layer index?
	#depth::Array{Float64, 1} # Depths for estimates 
	notes::String
end 

function Kohistan()
	deep = [(d > 52) & (d < 56) for d in crustDistribution.depth[:,4]]
	_, upper, middle, lower = mean(crustDistribution.depth[deep,:], dims=1) # layer cutoffs 

	## Garnet granulite gabbro starting at 41 km depth (from figure 4) + sarangar gabbro above 
	lower_comp = ((lower-41)/(lower-middle)) * 48.59 + (middle-41)/(middle-lower) * 52.52  

	## Remainder of SPC, normalized by geobarametric depth, minus the bit of the sarangar gabbro included in the lower crust 
	middle_spc_comp = ((4.3- (middle-41))*52.52 + 3.0*51.09 + .8*57.95 + 2.3*50.83 + .4*72.01 + 1.1*49.61) / ((4.3- (middle-41)) + 3.0 + .8 + 2.3 + .4 + 1.1)
	pluton_comp = (.07*53.64 + .67*64.49 + .09*73.37 ) / (.07 + .67 + .09)
	## SPC comp from 30 km (bottom of Gilgit complex, from table 1) to bottom of middle + pluton from upper boundary to 30km depth
	middle_comp = ((middle-30)*middle_spc_comp + (30-upper)*pluton_comp) / ((middle-30)+(30-upper))

	volcanics_comp = 56.83
	# 4km of volcanics (from table 1) + remainder pluton 
	upper_comp = (4*volcanics_comp + (upper-4)*pluton_comp)/upper

end 


results = [
	ResultWithError("This paper", [58.5, 64.8, 56.1], [.0317, .0298, .0223], [1,2,3], ""),
	Result("Rudnick and Fountain", [52.3, 60.6, 66.0], [3,2,1], "upper crust from Taylor and McLennan 1985"),
	Result("Rudnick and Gao", [53.4, 63.5, 66.6], [3,2,1], ""),
	Result("Hacker et al. (most felsic)", [62, 68], [3,2],
		"Most mafic felsic model of the three-layer models in this paper."),
	Result("Hacker et al. (most mafic)", [49, 53], [3,2],
		"Most mafic model of the three-layer models in this paper.")
]

results_comparerf = [
	ResultWithError("This paper", [58.5, 64.8, 56.1], [.0317, .0298, .0223], [1,2,3], ""),
	Result("Rudnick and Fountain", [52.3, 60.6, 66.0], [3,2,1], "upper crust from Taylor and McLennan 1985"),
	ResultWithError("This model with R&F data", [56.2,61.3,61.0], [.0468,.0399,.0422], [3,2,1], 
		"From visualization/compare_rf.jl run on latlong_nobin perplex dataset."),
	ResultWithError("This model with Crust1.0 data", [56.8, 62.4, 62.0], [.0382, .0434,.0439], [3,2,1], 
		"From visualization/compare_rf.jl run on latlong_nobin perplex dataset.")
]

# Updated as of 5/21/2020
results_all_versions = [
	ResultWithError("Base model", [58.5, 64.8, 56.1], [0.0303 , 0.0332 , 0.0336], [1,2,3], ""),
	ResultWithError("All earthchem samples, modified H&P", [54.9, 62.5, 53.8], [0.0279 , 0.0379 , 0.0329], [1,2,3], "remote/bulk_binned"),
	ResultWithError("Area-averaged Igneous samples", [59.7, 65.0, 58.1], [0.0289 , 0.0247 , 0.03], [1,2,3], "areaAve"),
	ResultWithError("All earthchem samples, base H&P", [55.3, 64.0, 56.3], [0.0276 , 0.037 , 0.0362], [1,2,3], "remote/bulk_binned_hp02"),
	ResultWithError("No geotherm binning", [59.6, 63.3, 58.6], [0.0383 , 0.0386 , 0.0405], [1,2,3], "remote/latlong_nobin"),
	ResultWithError("Igneous earthchem samples resampled to remove Daly gap", [58.4, 64.7, 56.6], [0.0267 , 0.0333 , 0.0342], [1,2,3],
		"remote/binned_geotherm_si_weighted"),
	#ResultWithError("TC1 age model", [58.5, 64.8, 56.1], [0.0303 , 0.0332 , 0.0336], [1,2,3], "remote/latlong_weighted, -a tc1"), # Exactly the same! Age model doesn't affect result.
	ResultWithError("3% seismic bin", [58.2, 64.1, 56.2], [0.0358 , 0.0321 , 0.0318], [1,2,3], "remote/latlong_weighted_bin3"),
	ResultWithError("10% seismic bin", [57.3, 63.7, 57.1], [0.031 , 0.0385 , 0.04], [1,2,3], "remote/latlong_weighted_bin10")
]

p = plot(xlabel="% SiO2", yflip=true, size=(300,500), legend=:outerbottom, legendtitlefontsize=9,
	fg_legend = :transparent, framestyle=:box)

for (r, result) in enumerate(results)
	if typeof(result) == Result
		plot!(result.composition, result.layers, label=result.name, markershape=:circle, markerstrokecolor=:auto)
	elseif typeof(result) == ResultWithError
		if r == 1
			plot!(result.composition, result.layers, xerror=(result.sem .* 2), 
				#color=:blue,linewidth=2,
				label=result.name, markershape=:circle, markerstrokecolor=:auto)
		else 
			plot!(result.composition, result.layers, xerror=(result.sem .* 2), 
				#color=:blue, markersize=2, linealpha=.5,
				label=result.name, markershape=:circle, markerstrokecolor=:auto)
		end
	else 
		println("Unrecognized result type")
	end 
end 

plot!(yticks=(1:3, ["Upper\ncrust","Middle\ncrust","Lower\ncrust"]))

savefig("output/comparison.pdf")




