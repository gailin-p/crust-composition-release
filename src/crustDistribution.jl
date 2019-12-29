"""
Load crust1 and tc1 data (layer depths and geotherms).
Weight each 1x1 degree by the actual size of that degree

"""
module crustDistribution

using StatsBase
using StatGeochem
using JLD

dataPath = "data/crust1Layers.jld"
lat_range = -89:90
long_range = -180:179
nlat = length(lat_range)
nlong = length(long_range)


depth = Array{Float64, 2} # n rows, 4 columns (upper, middle, lower, geotherm)
export depth
weights = Array{Float64, 1} # n rows, weight for choosing each row
export weights

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
    if isfile(dataPath)
        d = load(dataPath)
        global depth = d["depth"]
        global weights = d["weights"]
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
    contTest = cont .!= "NA"

    # Get depths.
    global depth = find_tc1_crust(lats, longs)
    for layer in [6,7,8] # upper, middle, lower
        d = [-1*n for n in find_crust1_base(lats, longs, layer)]
        depth = hcat(depth,d)
    end

    # Filter for only locations where all depths and geotherm are non-null
    test = [!isnan(x) for x in sum(depth,dims=2)] .& contTest
    depth = depth[test[:],:]
    lats = lats[test[:]] # latitudes for values we're keeping

    global weights = latLongWeight.(lats)

    save(dataPath,"depth",depth, "weights", weights)
end

"""
    getCrustParams()
Provide random layer depths drawn from distributions.
Return depths
    (550 isotherm, upper,middle,lower)
"""
function getCrustParams()
    i = sample(Weights(weights))
    return depth[i,:] # geotherm, base of upper, middle, lower
end
export getCrustParams

end
