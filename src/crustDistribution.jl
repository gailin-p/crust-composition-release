"""
Load crust1 and tc1 data (layer depths and geotherms).
Weight each 1x1 degree by the actual size of that degree

Usage note: run get_crust1() (in StatGeochem) to access Crust1.0 data.

TODO return dist values (gaussian, use uncertainty as std?) not exact

"""
module crustDistribution

using StatsBase
using StatGeochem
using JLD
using MAT
using Statistics
using Random

include("config.jl")
include("utilities.jl")
using Plots; gr();

# For resampling
#  If running from jupyter notebook, need relative path
uncertainty_dat = matread(IGN_FILE)["err2srel"]
relerr_tc1 = uncertainty_dat["tc1Crust"]/2

# 1 deg lat/long to account for uncertainty in placement of Crust1 boundaries.
err_latlong = 1.0

# Exports
depth = Array{Float64, 2} # n rows, 4 columns (geotherm, upper, middle, lower)
export depth
weights = Array{Float64, 1} # n rows, weight for choosing each row
export weights

# Lats and longs where both tc1 (so age) and crust1 are defined
all_lats = Array{Float64, 1} # n rows, age of lat/long for this row
all_longs = Array{Float64, 1} # n rows, age of lat/long for this row
export all_longs
export all_lats

# Ages of crust crust according to TC1
ages = Array{Float64, 1} # n rows
export ages

"""
    latLongWeight(lat, long)
return weighting for the area of a 1 degree square at this latitude.
"""
function latLongWeight(lat)
    lat = abs(lat)
    longLength = cos(deg2rad(lat))*69
    area = longLength * 69
    maxArea = 69^2
    return area/maxArea
end

"""
Load presaved dataset or write dataset
"""
function __init__()
    # Allow different relative paths so works in jupyter notebook or script.
    file = isdir("../resources") ? "../resources/crustDistribution.jld" : "resources/crustDistribution.jld"
    if isfile(file)
        loadSavedCrust(file)
    else
        loadCrust1()
        saveCrust(file)
    end
end

"""
    latsLongs()
Return lists of lats and longs for all 1x1 degree squares that are on continents.
Not all these values have both TC1 and Crust1 values, so will be futher filtered in loadCrust1.
"""
function latsLongs()
    # Valid range
    lat_range = -89:90
    long_range = -180:179
    nlat = length(lat_range)
    nlong = length(long_range)

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
    contTest = cont .!= "NA"
    longs = longs[contTest]
    lats = lats[contTest]

    ## Check Crust1.0. upper crust below 5.25 is ocean crust not caught by geocont
    (vp, vs, rho) = find_crust1_seismic(lats, longs, 6)
    vpTest = vp .> 5.25
    longs = longs[vpTest]
    lats = lats[vpTest]

    return (lats, longs)
end

function saveCrust(path)
    save(path, "weights", weights, "all_lats", all_lats, "all_longs", all_longs,
        "ages", ages, "depth", depth, "err_latlong", err_latlong)
end

function loadSavedCrust(path)
    println("Loading crust information from file $path")
    d = load(path)
    global weights = d["weights"]
    global all_lats = d["all_lats"]
    global all_longs = d["all_longs"]
    global ages = d["ages"]
    global depth = d["depth"]
    global err_latlong = d["err_latlong"]
end

"""Get thickness and depths
This doesn't work on Discovery (requires ImageCore, ImageMagick, which are broken)
So save to a file, then use __init__() to read that file
TODO: i think this is now fixed, so can dispense w this overhead
"""
function loadCrust1()
    (lats, longs) = latsLongs()
    #heatmap(globe(lats, longs, find_geolcont(lats, longs)))
    #savefig("before_filter.png")

    # Get depths.
    depth = find_tc1_crust(lats, longs) # Geotherm
    for layer in [6,7,8] # upper, middle, lower
        d = [-1*n for n in find_crust1_base(lats, longs, layer)]
        depth = hcat(depth,d)
    end

    # Filter for only locations where all depths and geotherm are non-null
    test = [!isnan(x) for x in sum(depth,dims=2)]
    depth = depth[test[:],:]
    lats = lats[test[:]] # latitudes for values we're keeping
    longs = longs[test[:]]
    ## Filtering by TC1 availibility signifcantly reduces the amount of
    #  continental shelf included in the lat/long pairs. 
    #heatmap(globe(lats, longs, find_geolcont(lats, longs)))
    #savefig("after_filter.png")

    global weights = latLongWeight.(lats)

    global all_lats = lats
    global all_longs = longs

    # Get ages.
    global ages = find_tc1_age(lats, longs)[1]

    global depth = depth

    global err_latlong = err_latlong
end


"""
    getCrustParams()
Provide random layer depths drawn from distributions.
Uses all geotherms (no geotherm binning)
Return depths
    (550 isotherm, upper,middle,lower)
"""
function getCrustParams(n::Int; uncertain::Bool=false)
    bottom = minimum(depth[:,1]) -1
    top = maximum(depth[:,1]) +1
    return getCrustParams(bottom, top, n, uncertain=uncertain)
end

# """
#     getFormationConditions(n)
# Return n depth to 550, depth pairs representing possible formation conditions.
# Units Bar, C
# Assume formation at base of crust.

# NOTE: After discussion with brenhin, this is not a good way to estimate formation conditions.
# Plutons do not form at base of crust; instead they generally form 0-10 km below where they are now.
# They form at ~600 C but experience some retrograding
# """
# function getFormationParams(n::Int)
#     test = ages .== 25 # youngest crust
#     samples = Array{Float64, 2}(undef, (n,2))
#     for i in 1:n
#         rand_i = sample(Weights(weights[test]))
#         this_sample = depth[test,:][rand_i,:] # depth to 550, base of upper, middle, lower
#         samples[i,1] = this_sample[1]
#         samples[i,2] = this_sample[4]
#         # # pressue is dpdz * base of lower crust
#         # samples[i,1] = dpdz * this_sample[4]
#         # # temperature: assume 0 at surface to 550 at isotherm
#         # dtdz = (550 / this_sample[1]) # deg c / km
#         # samples[i,2] = dtdz * this_sample[4]
#     end
#     return samples
# end

"""
    getFormationConditions()
Return n depth to 550, depth pairs representing mean formation conditions
(middle of crust at newly cratonized crust)

NOTE: After discussion with brenhin, this is not a good way to estimate formation conditions.
Plutons do not form at base of crust; instead they generally form 0-10 km below where they are now.
They form at ~600 C but experience some retrograding
"""
function getFormationParams()
    #test = ages .== 25 # youngest crust
    #depths = depth[test,:] ./ 2
    #w = Weights(weights[test])
    #return (mean(depths[:,1], w), mean(depths[:,4], w))
    t = 550 # degrees C
    d = 30 # km
    return (d, t/d)
end

# Return depth and dtdz of formation
function getFormationParams(depth::Number, dtdz::Number, uncertain::Bool=true)
    if uncertain
        fdepth = depth + rand()*10 # formation depth up to 10 km deeper than present depth
        ftemp = fdepth*dtdz + rand()*(600 - fdepth*dtdz) # formation temp between current geotherm and 600 C
    else
        fdepth = depth+10
        ftemp = 400
    end
    return fdepth, ftemp/fdepth
end

"""
    getCrustParams(min, max, n)
Get random layer depths and geotherms for n samples when geotherms are binned.
Use bin min and bin max to limit where geotherms are pulled from.
Bin prior to adding uncertainty, because that uncertainty in geotherm
exists for the samples we're inverting too.
Return Array of size (n, 4), geotherm and layers for n samples.
Used by resampleEarthChem to assign geotherms to resampled values.
"""
function getCrustParams(bin_min::Number, bin_max::Number, n::Int; uncertain::Bool=false)
    test = (depth[:,1] .<= bin_max) .& (depth[:,1] .> bin_min)
    bin_depth = depth[test,:]
    #println("bin has $(length(unique(bin_depth[:,1]))) samples") # Note that 88 is not in tc1 geotherms, so that bin will be smaller
    bin_weight = weights[test]
    samples = Array{Float64, 2}(undef, (n,4))
    for i in 1:n
        rand_i = sample(Weights(bin_weight))
        this_sample = bin_depth[rand_i,:] # geotherm, base of upper, middle, lower
        if uncertain
            tc1Std = this_sample[1]*relerr_tc1 # standard deviation is relative to geotherm
            crustStd = this_sample[3]*relerr_tc1 # Use middle crust to get std TODO ???
            samples[i,:] = [this_sample[1]+(randn()*tc1Std), this_sample[2]+(randn()*crustStd),
                this_sample[3]+(randn()*crustStd), this_sample[4]+(randn()*crustStd)]
        else
            samples[i,:] = this_sample
        end
    end
    return samples
end
export getCrustParams

"""
    binBoundaries(n)
Return range of boundaries that gives x bins
Ensure that range is divisible by n by padding on top and bottom
so that most bins have same # of discrete tc1 geotherm depths,
so that resulting geotherm ranges in data sent to perplex will be of equal size
"""
function binBoundaries(n::Int)
    dmin = minimum(depth[:,1])
    dmax = maximum(depth[:,1])
    return range(dmin, dmax, length=n+1)
end
export binBoundaries

end
