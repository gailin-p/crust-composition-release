"""
    commutativityOfAverage
Functions for selecting random batches of samples and running perplex on them.
"""
module commutativityOfAverage

using ProgressMeter
using StatGeochem
using Statistics
using DelimitedFiles
using Random
using JLD

include("crustDistribution.jl")

# Average geotherm and layer depths
# General options
exclude = ""
dataset = "hpha02ver.dat"
elts = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]
dpdz = 2900. * 9.8 / 1E5 * 1E3
solutions_nof = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\n"
npoints = 20
# upperStart = 1.08962179685088
# upperBase = 13.442666718122878
# middleBase = 25.551929607903677
# lowerBase = 36.57912627354122
# depth = 40.2 # mean 550C depth of all samples

## Set up variables
# Using Crust1.0 interface provided by StatGeochem, get all layer depths
layers = Dict{String,Int}(("Upper"=>6),("Middle"=>7),("Lower"=>8)) # Crust1.0 layers

"""
    badValToNan!(a)
Convert values in a < lower or > upper to NaN
"""
function badValToNan!(a::Array{Float64,1}; lower::Float64=1e-6, upper::Float64=1e6)
    # Ignore very small and very large values (we'll use nanmean)
    badval = map(s-> (s > upper || s < lower),a)
    if (sum(badval)) > 0
        println("$(sum(badval)) bad values")
    end
    a[badval] .= NaN
    return a
end

# TODO: This is the same as massMean, is that correct?
function seismicMean(nums)
    oneover = [1/n for n in nums]
    return 1/nanmean(oneover)
end
export seismicMean

function massMean(nums)
    # n = length(densities)
    # volumes = [(1/d)*(1/n) for d in densities]
    # return 1/sum(volumes)
    oneover = [1/n for n in nums]
    return 1/nanmean(oneover)
end
export massMean

"""
    runSamples!()
Run r samples
Return (ave properties, properties of averages, indices)
properties arrays are 3x3xr, layer x property x run
TODO refactor to return new instead of modify result arrays
"""
function runSamples!(props_of_ave::Array{Float64,3}, ave_properties::Array{Float64,3}, indices::Array{Int64,2},
    r::Int, n::Int, ignFile::String, perplex::String, scratch::String)
    dat = readdlm(ignFile,',')

    myindex = rand(Int)

    # run r sample pairs through perplex
    @showprogress 1 "Running samples... " for i in 1:r
        # For each sample pair, go get a new crust depth_hist
        (tc1, upperBase, middleBase, lowerBase) = crustDistribution.getCrustParams()
        geotherm = 550.0/tc1/dpdz
        P_range = [1, ceil(Int,lowerBase*dpdz)] # run to base of crust. TODO should we only run perplex to 550 degree isotherm?

        allComp = Array{Float64,2}(undef,(n,10))
        allSeismic = Array{Any,1}(undef,n)

        # Run perplex for each composition
        for j in 1:n
            index = convert(Int, ceil(rand()*size(dat)[1]))
            indices[j,i] = index
            sample = dat[index,:]
            comp = sample[2:11]
            allComp[j,:] = comp
            # Run perple_x
            perplex_configure_geotherm(perplex, scratch, comp, elements=elts,
                P_range=P_range, geotherm=geotherm, dataset=dataset, solution_phases=solutions_nof,
                excludes="", index=myindex, npoints=npoints)

            seismic = perplex_query_seismic(perplex, scratch, index=myindex)
            # Ignore very small and very large values
            badValToNan!(seismic["rho,kg/m3"])
            badValToNan!(seismic["vp,km/s"],lower=1.0)
            badValToNan!(seismic["vp/vs"])
            allSeismic[j]=seismic
        end

        # Average seismic properties over samples (assume p is same at each index)
        # keys: "vp,km/s", "rho,kg/m3", "T(K)", "elements", "P(bar)", "vp/vs"
        aveVp = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
        aveVpVs = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
        aveRho = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))

        for p_index in 1:length(allSeismic[1]["P(bar)"])
            aveVp[p_index] = seismicMean(map(y -> (allSeismic[y]["vp,km/s"][p_index]), 1:length(allSeismic)))
            aveVpVs[p_index] = seismicMean(map(y -> (allSeismic[y]["vp/vs"][p_index]), 1:length(allSeismic)))
            aveRho[p_index] = massMean(map(y -> (allSeismic[y]["rho,kg/m3"][p_index]), 1:length(allSeismic)))
        end

        # Run perplex on average composition
        aveComp = Array{Float64,1}(undef,10)
        for j in 1:10
            aveComp[j] = mean(allComp[:,j])
        end
        perplex_configure_geotherm(perplex, scratch, aveComp, elements=elts,
            P_range=P_range, geotherm=geotherm, dataset=dataset, solution_phases=solutions_nof,
            excludes="", npoints=npoints, index=myindex)

        seismicOfAve = perplex_query_seismic(perplex, scratch, index=myindex)
        # Ignore very small and very large values
        badValToNan!(seismicOfAve["rho,kg/m3"])
        badValToNan!(seismicOfAve["vp,km/s"],lower=4.0) # See if this clears up the bimodality
        badValToNan!(seismicOfAve["vp/vs"])

        # Average seismic properties over depth layer
        # (TODO is this ok? we're trying to *test* the legitimacy of averaging)
        upperTest = seismicOfAve["P(bar)"] .< upperBase * dpdz
        middleTest = (seismicOfAve["P(bar)"] .> upperBase * dpdz) .& (seismicOfAve["P(bar)"] .< middleBase * dpdz)
        lowerTest = (seismicOfAve["P(bar)"] .> middleBase * dpdz) .& (seismicOfAve["P(bar)"] .< lowerBase * dpdz)
        tests = [upperTest, middleTest, lowerTest]

        for l in 1:3 # upper, middle, lower
            ave_properties[l,1,i] = nanmean(aveRho[tests[l]])
            ave_properties[l,2,i] = nanmean(aveVp[tests[l]])
            ave_properties[l,3,i] = nanmean(aveVpVs[tests[l]])
            props_of_ave[l,1,i] = nanmean(seismicOfAve["rho,kg/m3"][tests[l]])
            props_of_ave[l,2,i] = nanmean(seismicOfAve["vp,km/s"][tests[l]])
            props_of_ave[l,3,i] = nanmean(seismicOfAve["vp/vs"][tests[l]])
        end
    end

    return (props_of_ave, ave_properties, indices)
end
export runSamples!

end
