"""
Using three group sizes created by commutativityOfAveragingMPI.jl, create histograms of difference for each. 
"""

using Plots; gr(); 
using JLD

files = ["data/commutativity-n2-r2007.jld", "data/commutativity-n10-r2007.jld", "data/commutativity-n100-r2007.jld"]

plots = []

lims = (-.15, .15)

for (i, prop) in enumerate(["Rho", "Vp", "Vp/Vs"])
	d = load(files[1])
	diffs = (d["props_of_ave"][:,i,:] .- d["ave_properties"][:,i,:]) ./ d["ave_properties"][:,i,:]
	p = plot(collect(Iterators.flatten(diffs)), normalize=:pdf, seriestype=:stephist, label="N=2", 
		xlims=lims, yaxis=false, grid=false)

	d = load(files[2])
	diffs = (d["props_of_ave"][:,i,:] .- d["ave_properties"][:,i,:]) ./ d["ave_properties"][:,i,:]
	p = plot!(collect(Iterators.flatten(diffs)), normalize=:pdf, seriestype=:stephist, label="N=10", 
		xlims=lims, yaxis=false, grid=false)

	d = load(files[3])
	diffs = (d["props_of_ave"][:,i,:] .- d["ave_properties"][:,i,:]) ./ d["ave_properties"][:,i,:]
	p = plot!(collect(Iterators.flatten(diffs)), normalize=:pdf, seriestype=:stephist, label="N=100", 
		xlims=lims, yaxis=false, grid=false)

	p = plot!(title=prop)

	if i != 1
		p = plot!(legend=false)
	end 

	append!(plots, [p])
end 

layer_plots = []
for (i, layer) in enumerate(["Upper", "Middle", "Lower"])
	d = load(files[1])
	diffs = (d["props_of_ave"][i,2,:] .- d["ave_properties"][i,2,:]) ./ d["ave_properties"][i,2,:]
	p = plot(collect(Iterators.flatten(diffs)), normalize=:pdf, seriestype=:stephist, label="N=2", 
		xlims=lims, yaxis=false, fg_legend = :transparent, grid=false)

	d = load(files[2])
	diffs = (d["props_of_ave"][i,2,:] .- d["ave_properties"][i,2,:]) ./ d["ave_properties"][i,2,:]
	p = plot!(collect(Iterators.flatten(diffs)), normalize=:pdf, seriestype=:stephist, label="N=10", 
		xlims=lims, yaxis=false, grid=false)

	d = load(files[3])
	diffs = (d["props_of_ave"][i,2,:] .- d["ave_properties"][i,2,:]) ./ d["ave_properties"][i,2,:]
	p = plot!(collect(Iterators.flatten(diffs)), normalize=:pdf, seriestype=:stephist, label="N=100", 
		xlims=lims, yaxis=false, grid=false)

	p = plot!(title=layer)

	if i != 1
		p = plot!(legend=false)
	end

	if i == 3
		p = plot!(xlabel="Difference as % of average Vp")
	end

	append!(layer_plots, [p])
end 

p = plot(plots..., layout=(3,1), size=(400, 500))
savefig(p, "output/commutativity_hist_by_prop.pdf")

p2 = plot(layer_plots..., layout=(3,1), size=(400, 500))
savefig(p2, "output/commutativity_hist_by_layer.pdf")



