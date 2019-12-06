module commutativityOfAverage

using ProgressMeter
using StatGeochem
using Statistics
using DelimitedFiles
using Random
using JLD

# Assumes Crust1.0 data as produced by loadAndSaveCrust1():
dataPath = "data/crust1Layers.jld"

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
lat_range = -89:90
long_range = -180:179
nlat = length(lat_range)
nlong = length(long_range)
layers = Dict{String,Int}(("Upper"=>6),("Middle"=>7),("Lower"=>8)) # Crust1.0 layers
depth = Array{Float64, 2} # n rows, 4 columns (upper, middle, lower, geotherm)
export depth

"""
Load and save doesn't work on discovery
"""
function __init__()
    if isfile(dataPath)
        d = load(dataPath)
        global depth = d["depth"]
    else
        loadAndSaveCrust1()
    end
end

"""Get thickness and depths
This doesn't work on Discovery (requires ImageCore, ImageMagick, which are broken)
So save to a file, then use __init__() to read that file.
"""
function loadAndSaveCrust1()
    longs = Array{Float64,1}()
    lats = Array{Float64,1}()
    for lat in lat_range
        for long in long_range
            push!(longs, long)
            push!(lats, lat)
        end
    end

    ## Check continent
    cont = map(c -> continents[c], find_geolcont(lats,longs))
    test = cont .== "NA"
    depthTest = [] # for filtering isotherm by depth NaNs

    # Get depths.
    global depth = find_tc1_crust(lats, longs)
    for layer in [6,7,8] # upper, middle, lower
        d = [-1*n for n in find_crust1_base(lats, longs, layer)]
        depth = hcat(depth,d)
    end

    # Filter for only locations where all depths and geotherm are non-null
    test = [!isnan(x) for x in sum(depth,dims=2)]
    depth = depth[test[:],:]

    save(dataPath,"depth",depth)
end

"""
    getCrustParams()
Provide random layer depths drawn from distributions.
Return depths
    (upper,middle,lower)
"""
function getCrustParams()
    i = ceil(Int,size(depth)[1]*rand())
    return depth[i,:] # geotherm, base of upper, middle, lower
end

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

"""
    runSamples!()
Run r samples
Return (ave properties, properties of averages, indices)
properties arrays are 3x3xr, layer x property x run
"""
function runSamples!(props_of_ave::Array{Float64,3}, ave_properties::Array{Float64,3}, indices::Array{Int64,2},
    r::Int, n::Int, ignFile::String, perplex::String, scratch::String)
    dat = readdlm(ignFile,',')

    myindex = rand(Int)

    # run r sample pairs through perplex
    @showprogress 1 "Running samples... " for i in 1:r
        # For each sample pair, go get a new crust depth_hist
        (tc1, upperBase, middleBase, lowerBase) = getCrustParams()
        print("$tc1, $upperBase, $middleBase, $lowerBase")
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
            badValToNan!(seismic["vp,km/s"])
            badValToNan!(seismic["vp/vs"])
            allSeismic[j]=seismic
        end

        # Average seismic properties over samples (assume p is same at each index)
        # keys: "vp,km/s", "rho,kg/m3", "T(K)", "elements", "P(bar)", "vp/vs"
        aveVp = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
        aveVpVs = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
        aveRho = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))

        for p_index in 1:length(allSeismic[1]["P(bar)"])
            aveVp[p_index] = nanmean(map(y -> (allSeismic[y]["vp,km/s"][p_index]), 1:length(allSeismic)))
            aveVpVs[p_index] = nanmean(map(y -> (allSeismic[y]["vp/vs"][p_index]), 1:length(allSeismic)))
            aveRho[p_index] = nanmean(map(y -> (allSeismic[y]["rho,kg/m3"][p_index]), 1:length(allSeismic)))
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
