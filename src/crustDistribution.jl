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
relerr_tc1 = uncertainty_dat["tc1Crust"]/2
relerr_crust = uncertainty_dat["Crust"]/2 # for crust1.0

# Exports 
depth = Array{Float64, 2} # n rows, 4 columns (upper, middle, lower, geotherm)
export depth
weights = Array{Float64, 1} # n rows, weight for choosing each row
export weights
ages = Array{Float64, 1} # n rows, age of lat/long for this row 
export ages

# Lats and longs where both tc1 (so age) and crust1 are defined 
all_lats = Array{Float64, 1} # n rows, age of lat/long for this row 
all_longs = Array{Float64, 1} # n rows, age of lat/long for this row 


all_ages = [ NaN  NaN  NaN
                25    0   50
               150   50  250
               395  250  540
               695  540  850
               975  850 1100
              1400 1100 1700
              2100 1700 2500
              2750 2500 3000
              3250 3000 3500] # From tc1.jl, bins for ages of crust according to tc1 

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
        global ages = d["ages"]
        global all_lats = d["lats"]
        global all_longs = d["longs"]
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
    longs = longs[test[:]] 

    # get ages 
    global ages = find_tc1_age(lats, longs)[1]

    global weights = latLongWeight.(lats)

    global all_lats = lats 
    global all_longs = longs 

    save(dataPath,"depth",depth, "weights", weights, "ages", ages, "lats", lats, "longs", longs)
end

"""
    getCrustParams()
Provide random layer depths drawn from distributions.
Return depths
    (550 isotherm, upper,middle,lower)
TODO return dist values (gaussian, use uncertainty as std?) not exact
"""
function getCrustParams(n::Int; uncertain::Bool=false)
    samples = Array{Float64, 2}(undef, (n,4))
    for i in 1:n
        rand_i = sample(Weights(weights))
        this_sample = depth[rand_i,:] # geotherm, base of upper, middle, lower
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


"""
    getCrustParams(bin, max, n)
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
    to_add = (n - 1 - (dmax - dmin)%(n-1))/2 # add this to top and bottom so that (dmax - dmin)%(n-1) = 0 and we have even-sized bins
    return LinRange(dmin-to_add, dmax+to_add, n)
end
export binBoundaries 

"""
    getAllSeismic()
Return a touple of weights and lists of seismic values (rho, vp, vp/vs)
From Crust1.0
Also resample (using weights corresponding to area of each 1x1 degree square).
"""
function getAllSeismic(layer::Integer; age::Number=NaN, n=50000)
    if !(layer in [6,7,8])
        throw(ArgumentError("Layer must be 6, 7, or 8 (crysteline Crust1 layers)"))
    end
    if !((age in all_ages[:,1]) | isnan(age))
        throw(ArgumentError("Age must be one of the ages returned by tc1"))
    end

    global all_lats # Lats where both tc1 and crust1.0 defined 
    global all_longs
    global depth
    global ages

    if !isnan(age) # filter for this age bin 
        test = ages .== age
        lats = all_lats[test]
        longs = all_longs[test]
        geotherms = depth[test,4]
    else 
        lats = all_lats 
        longs = all_longs 
        geotherms = depth[:,4]
    end 

    (vp, vs, rho) = find_crust1_seismic(lats, longs, layer)

    k = latLongWeight.(lats)
    # Probability of keeping a given data point when sampling:
    # We want to select roughly one-fith of the full dataset in each re-sample,
    # which means an average resampling probability <p> of about 0.2
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)

    # Resample! 
    samples = hcat(rho, vp, vs, geotherms)
    sigma = hcat(rho .* relerr_rho, vp .* relerr_vp, vs .* relerr_vs, geotherms .* relerr_tc1)
    resampled = bsresample(samples, sigma, n, p)

    rho = resampled[:,1]
    vp = resampled[:,2]
    vs = resampled[:,3]
    geotherms = resampled[:,4]

    return (rho, vp, vp ./ vs, geotherms)
end

end
