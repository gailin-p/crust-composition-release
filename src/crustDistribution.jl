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

# For building querries 
dataPath = "data/crust1Layers.jld"
lat_range = -89:90
long_range = -180:179
nlat = length(lat_range)
nlong = length(long_range)

# For resampling 
uncertainty_dat = matread("igncn1.mat")["err2srel"]
relerr_rho = uncertainty_dat["Rho"]/2
relerr_vp = uncertainty_dat["Vp"]/2
relerr_vs = uncertainty_dat["Vs"]/2

# Exports 
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

"""
    latsLongs() 
Return lists of lats and longs for all 1x1 degree squares that are on continents
"""
function latsLongs() 
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

    return (lats, longs)
end

"""Get thickness and depths
This doesn't work on Discovery (requires ImageCore, ImageMagick, which are broken)
So save to a file, then use __init__() to read that file.
"""
function loadAndSaveCrust1()
    (lats, longs) = latsLongs()

    # Get depths.
    global depth = find_tc1_crust(lats, longs)
    for layer in [6,7,8] # upper, middle, lower
        d = [-1*n for n in find_crust1_base(lats, longs, layer)]
        depth = hcat(depth,d)
    end

    # Filter for only locations where all depths and geotherm are non-null
    test = [!isnan(x) for x in sum(depth,dims=2)] 
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
TODO return dist values (gaussian, use uncertainty as std?) not exact
"""
function getCrustParams()
    i = sample(Weights(weights))
    return depth[i,:] # geotherm, base of upper, middle, lower
end
export getCrustParams

"""
    getAllSeismic()
Return a touple of weights and lists of seismic values (weights, rho, vp, vp/vs)
Weights are for differing sizes of lat/long squares 
  (different from global weights because different filter)
Also resample. 
"""
function getAllSeismic(layer::Int64)
    if !(layer in [6,7,8])
        throw(ArgumentError("Layer must be 6, 7, or 8 (crysteline Crust1 layers)"))
    end

    (lats, longs) = latsLongs()

    (vp, vs, rho) = find_crust1_seismic(lats, longs, layer)

    k = latLongWeight.(lats)
    # Probability of keeping a given data point when sampling:
    # We want to select roughly one-fith of the full dataset in each re-sample,
    # which means an average resampling probability <p> of about 0.2
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)

    # Resample! 
    samples = hcat(rho, vp, vs)
    sigma = hcat(rho .* relerr_rho, vp .* relerr_vp, vs .* relerr_vs)
    resampled = bsresample(samples, sigma, 50000, p)

    rho = resampled[:,1]
    vp = resampled[:,2]
    vs = resampled[:,3]

    return (rho, vp, vp ./ vs)
end
export getAllSeismic

end
