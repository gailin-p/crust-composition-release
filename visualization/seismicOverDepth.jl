## For one perplex output sample, plot seismic properties over depth

using Plots; gr();
using StatGeochem

I = "local" # index to use

dpdz = 2900 * 9.8 / 1E5 * 1E3 # Pressure gradient (bar/km). From runPerplexBatchVp.c

perplexresults = importdataset("../data/runTestNew.small.log", '\t')
#perplexresults = importdataset("../../brenhin_code/inversion/PerplexResults.NoMelt.log", '\t')

all = 1:1:length(perplexresults["P(bar)"])
test = (perplexresults["index"] .== I) .& ((all .% 10) .== 0) # index of sample to plot

depth = perplexresults["P(bar)"][test] ./ dpdz

p1 = plot(perplexresults["rho"][test], depth, xflip=true, xlabel="Rho", rotation=45, legend=false, ylabel="Depth (km)")
p2 = plot(perplexresults["Vp"][test], depth, xflip=true, xlabel="Vp", yaxis=nothing, legend=false)
p3 = plot(perplexresults["Vp/Vs"][test], depth, xflip=true, xlabel="Vp/Vs", yaxis=nothing, legend=false)

p = plot(p1, p2, p3, layout=(1,3), size=(800,800))

savefig(p, "../output/seismic-over-depth-"*string(I)*".pdf")
