### For some number of random indices of samples from ignmajors.csv,
# run perplex on the average composition of those samples and on those samples individually.
# Plot seismic properties over depth for averages vs average of seismic properties
# Just figure generation, not test of average effect

# Pairs that create bimodal distribution:
# [27065, 9246][21002, 41595][52903, 23546][32045, 666][5049, 25333][23325, 43253]
#[7983, 6440][14806, 6399][30548, 9059][30611, 42222][30886, 26707][47774, 27331][3085, 4350]
#[26289, 42503][40559, 19392][57807, 3327][28725, 31168][3320, 20499][5177, 17946]
#[46927, 33603][50768, 5046][45112, 3302][5453, 9252][32847, 13125][48947, 22786]
#[6440, 35041][15498, 34660][23323, 2798][45015, 14809][9227, 40488][22932, 6303][21173, 53182]
#[48191, 10666][48267, 8486][29917, 12136][34998, 22786][1723, 21167][9109, 13199][42326, 28655]
#[25768, 32071][56476, 29061][53752, 42398][44658, 41588][16985, 46258][26897, 21002][9497, 44658]
#[10036, 24390][31926, 42717][38113, 5194][17387, 10020]

using DelimitedFiles
using Random
using Statistics
using StatGeochem
using Plots; gr();
using ProgressMeter
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--indices"
        help = "Indices of samples to run"
        arg_type = Int
        nargs = '*'
        required = true
end
parsed_args = parse_args(ARGS, s)

# General options
perplexdir = "/Users/gailin/dartmouth/crustal_structure/Perple_X/"
scratchdir = "/Users/gailin/dartmouth/crustal_structure/perplexed_pasta/"
exclude = ""
dataset = "hpha02ver.dat"
elts = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]
P_range = [1, 20000]
dpdz = 2900. * 9.8 / 1E5 * 1E3
depth = 40.2 # mean 550C depth of all samples
geotherm = 550.0/depth/dpdz
solutions_nof = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\n"
solutions_f_melt = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\nF\nmelt(HP)\n"
npoints = 30

# Where to put figures
prefix = "../output/perplex/"

dat = readdlm("../ignmajors.csv",',')

count = length(parsed_args["indices"])

function seismicMean(nums)
    oneover = [1/n for n in nums]
    return 1/nanmean(oneover)
end

function badValToNan!(a::Array{Float64,1}; lower::Float64=1e-6, upper::Float64=1e6)
    # Ignore very small and very large values (we'll use nanmean)
    badval = map(s-> (s > upper || s < lower),a)
    if (sum(badval)) > 0
        println("$(sum(badval)) bad values")
    end
    a[badval] .= NaN
    return a
end

# RUN PERPLEX

allComp = Array{Float64,2}(undef,(count,10))
allSeismic = Array{Any,1}(undef,count)

# Run perplex for each composition
@showprogress 1 "Running samples: " for (i,index) in enumerate(parsed_args["indices"])
    sample = dat[index,:]
    comp = sample[2:11]
    allComp[i,:] = comp
    # Run perple_x
    perplex_configure_geotherm(perplexdir, scratchdir, comp, elements=elts,
        P_range=P_range, geotherm=geotherm, dataset=dataset, solution_phases=solutions_nof,
        excludes="h2oL\n", npoints=npoints)

    seismic = perplex_query_seismic(perplexdir, scratchdir)

    # Ignore very small and very large values
    badValToNan!(seismic["rho,kg/m3"])
    badValToNan!(seismic["vp,km/s"],lower=4.0)
    badValToNan!(seismic["vp/vs"])

    allSeismic[i]=seismic
end

# Average seismic properties
# keys: "vp,km/s", "rho,kg/m3", "T(K)", "elements", "P(bar)", "vp/vs"
for i in 1:count # Check that seismic profiles are equivalent
    if allSeismic[1]["T(K)"] != allSeismic[i]["T(K)"]
        println("calced geotherm for $(i) not the same")
        println(allSeismic[1]["T(K)"][1:10])
        println(allSeismic[i]["T(K)"][1:10])
    end
end

# Average seismic properties over samples (assume p is same at each index)
# keys: "vp,km/s", "rho,kg/m3", "T(K)", "elements", "P(bar)", "vp/vs"
aveVp = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
aveVpVs = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
aveRho = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))

for p_index in 1:length(allSeismic[1]["P(bar)"])
    aveVp[p_index] = seismicMean(map(y -> (allSeismic[y]["vp,km/s"][p_index]), 1:length(allSeismic)))
    aveVpVs[p_index] = seismicMean(map(y -> (allSeismic[y]["vp/vs"][p_index]), 1:length(allSeismic)))
    aveRho[p_index] = seismicMean(map(y -> (allSeismic[y]["rho,kg/m3"][p_index]), 1:length(allSeismic)))
end

# Plot averaged seismic properties
depth = allSeismic[1]["P(bar)"] ./ dpdz


# For now, just plot vp
p2 = plot(aveVp, size=(200,400), depth, yflip=true, xlabel="Vp (km/s)", ylabel="Depth (km)", label="ave vp", legend=:bottomleft);
for i in 1:length(allSeismic)
    plot!(p2, allSeismic[i]["vp,km/s"], depth, linestyle=:dot,
        label = "sample $(string(parsed_args["indices"][i]))")
    println("sample $(string(parsed_args["indices"][i])) has minimum vp $(nanminimum(allSeismic[i]["vp,km/s"]))")
    nans = [(x<1e-6) for x in allSeismic[i]["vp,km/s"]]
    println("sample $(string(parsed_args["indices"][i])) has $(sum(nans)) zeros out of $(length(allSeismic[i]["vp,km/s"]))")
end

#p1 = plot(aveRho, depth, xflip=true, xlabel="Rho (kg/m3)", color=:blue, rotation=45, legend=true, ylabel="Depth (km)", label="ave properties");
#plot!(p1, allSeismic[1]["rho,kg/m3"], depth, linestyle=:dot, color=:blue, alpha=0.5, label="first comp")
#plot!(p1, allSeismic[2]["rho,kg/m3"], depth, linestyle=:dot, color=:blue, alpha=0.5, label="second comp")
#p2 = plot(aveVp, depth, xflip=true, xlabel="Vp (km/s)", color=:blue, yaxis=nothing, legend=false);
#plot!(p2, allSeismic[1]["vp,km/s"], depth, linestyle=:dot, color=:blue, alpha=0.5)
#plot!(p2, allSeismic[2]["vp,km/s"], depth, linestyle=:dot, color=:blue, alpha=0.5)
#p3 = plot(aveVpVs, depth, xflip=true, xlabel="Vp/Vs", color=:blue, yaxis=nothing, legend=false);
#plot!(p3, allSeismic[1]["vp/vs"], depth, linestyle=:dot, color=:blue, alpha=0.5)
#plot!(p3, allSeismic[2]["vp/vs"], depth, linestyle=:dot, color=:blue, alpha=0.5)

# Run perplex on average composition
aveComp = Array{Float64,1}(undef,10) # number of comp elements
for i in 1:10
    aveComp[i] = mean(allComp[:,i])
end
perplex_configure_geotherm(perplexdir, scratchdir, aveComp, elements=elts,
    P_range=P_range, geotherm=geotherm, dataset=dataset, solution_phases=solutions_nof,
    excludes="", npoints=npoints)

seismicOfAve = perplex_query_seismic(perplexdir, scratchdir)

# # Plot seismic properties of average
#plot!(p1,seismicOfAve["rho,kg/m3"], depth,label="properties of ave")
plot!(p2,seismicOfAve["vp,km/s"], depth, label="average comp")
#plot!(p3,seismicOfAve["vp/vs"], depth)

plot!(p2, [1,1,1], [13.44, 25.55, 36.58], seriestype=:scatter, marker=7,markershape=:star5, label="crust layer boundaries")

#p = plot(p1, p2, p3, layout=(1,3), size=(800,600));
#savefig(p, "$(prefix)average-seismic-properties-$(count)samples$(suffix).pdf");
label = join([string(n) for n in sort(parsed_args["indices"])], "-")
savefig(p2, "$(prefix)average-seismic-properties-$(label).pdf")
