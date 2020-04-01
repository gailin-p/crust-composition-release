"""
Data from Crust1.0 to invert. 

Allows choice of age models. 
"""

include("crustDistribution.jl")
include("config.jl")

using MAT
using StatGeochem
using StatsBase
using Statistics
using Test 
using Random

uncertainty_dat = matread("igncn1.mat")["err2srel"]
relerr_rho = uncertainty_dat["Rho"]/2
relerr_vp = uncertainty_dat["Vp"]/2
relerr_vs = uncertainty_dat["Vs"]/2
relerr_tc1 = uncertainty_dat["tc1Crust"]/2
relerr_crust = uncertainty_dat["Crust"]/2 # for crust1.0

# TODO what's a reasonable lat/long error for crust 1.0? For now, use .5 deg lat/long 
err_latlong = .5

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


abstract type AgeModel end 

struct Tc1Age <: AgeModel end 

mutable struct EarthChemAge <: AgeModel 
	nbins::Integer # Number of age bins 
	margin::Integer # How many (equator-sized) blocks (outside of this one) to look in for an age 

	# For lookup 
	lats::Array{Float64, 1}
	longs::Array{Float64, 1}
	ages::Array{Float64, 1}
end 

function EarthChemAge(nbins::Integer, margin::Integer)
	println("Using EarthChem age model")
	ign = matread(IGN_FILE)
	return EarthChemAge(nbins, margin, ign["Latitude"][:], ign["Longitude"][:], ign["Age"][:])
end 

function ageBins(ageModel::EarthChemAge)
	return range(0, length=ageModel.nbins+1, stop=4000)
end 

function sampleAge(lats::Array{Float64, 1}, longs::Array{Float64, 1}, model::EarthChemAge)
	means = fill(-1.0, length(lats))
	uncertainties = fill(-1.0, length(lats))
	for (l, lat) in enumerate(lats)
		longmargin = model.margin / crustDistribution.latLongWeight(lat) # weight for this lat: larger margin in lat/long as size decreases
       	long = longs[l] 
       	test = (abs.(model.lats .- lat) .< model.margin) .& (abs.(model.longs .- long) .< longmargin) # latitude stays the same 
       	if sum(test) >= 1
       		index = ceil(Int, rand()*sum(test))
       		means[l] = model.ages[test][index] 
       		uncertainties[l] = nanstd(model.ages[test])
       	#	means[l] = maximum(model.ages[test])
       	#	means[l] = nanmean(model.ages[test])
       	end
    end
    #@test minimum(means) > 0 
    # TODO nan values should use TC1 age 
    bad = isnan.(means) .| (means .< 0)
    println("Using $(sum(bad)) TC1 values for ages")
    means[bad] .= find_tc1_age(lats[bad], longs[bad])[1]

    # Fix bad uncertainties (may be more than bad ages for locations with only one sample)
    bad = isnan.(uncertainties) .| (uncertainties .< 0)
    tcages = find_tc1_age(lats[bad], longs[bad])
    uncertainties[bad] .= (tcages[3] .- tcages[2]) ./ 2 

    return means, uncertainties
end 

""" TC1 age model """

# Return enumerable list of age bins (tops and bottoms)
function ageBins(ageModel::Tc1Age)
	return vcat(all_ages[2:end,2], [all_ages[end,3]]) # Bins 
end 

function sampleAge(lats::Array{Float64, 1}, longs::Array{Float64, 1}, ageModel::Tc1Age)
	sampleAges = find_tc1_age(lats, longs) # center, bin min, bin max 
	sampleUncertainties = (sampleAges[3] .- sampleAges[2]) ./ 2 

	if minimum(sampleUncertainties) < 0 
		throw(AssertionError("Uncertainty less than 0"))
	end 

	return sampleAges[1], sampleUncertainties
end 

"""
  nOriginal 
Get number of seismic data points prior to resample 
"""
function nOriginal()
  return length(crustDistribution.all_lats)
end 

"""
  FRSeismic() 
Same idea as getAllSeismic, but return fountain/rudnick data. Only returns vp. 
layer: 6, 7, or 8 (upper, middle, lower)
"""
function FRSeismic(layer::Integer; model::Integer=2, n::Integer=400, resample::Bool=true)
  dat, header = readdlm("Fountain-Rudnick.csv", ',', header=true)
  header = header[:]

  # Get data 
  layers = ["upper vp", "middle vp", "lower vp"]
  vp = dat[:,findfirst(isequal(layers[layer-5]), header)]
  weights = dat[:,findfirst(isequal("%, RF model $(model)"), header)]
  good = vp .!= ""
  vp = float.(vp[good])
  weights = float.(weights[good])

  # Resample 
  if resample 
    resampled = bsresample(vp, relerr_vp .* vp, n, weights)
    return resampled[:] 
  else 
    return vp 
  end 
end 


"""
    getAllSeismic()
Return a touple of weights and lists of seismic values (rho, vp, vp/vs, geotherm)
From Crust1.0
Optionally resample (using weights corresponding to area of each 1x1 degree square).
Latlong info optionally returned. 
By default, resample with Tc1 ages 
"""
function getAllSeismic(layer::Integer; ageModel::AgeModel=Tc1Age(), n::Integer=50000, resample::Bool=true, latlong::Bool=false)
    if !(layer in [6,7,8])
        throw(ArgumentError("Layer must be 6, 7, or 8 (crysteline Crust1 layers)"))
    end

    # Gather data... 
    lats = crustDistribution.all_lats # Lats where both tc1 and crust1.0 defined 
    longs = crustDistribution.all_longs
    geotherms = crustDistribution.depth[:,1]
    crustbase = crustDistribution.depth[:,4]

    ages, age_uncertainty = sampleAge(lats, longs, ageModel)

    (vp, vs, rho) = find_crust1_seismic(lats, longs, layer)

    if resample 
        k = crustDistribution.latLongWeight.(lats)
        # Probability of keeping a given data point when sampling:
        # We want to select roughly one-fith of the full dataset in each re-sample,
        # which means an average resampling probability <p> of about 0.2
        p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)

        if latlong 
            samples = hcat(rho, vp, vs, geotherms, crustbase, ages, lats, longs)
            sigma = hcat(rho .* relerr_rho, vp .* relerr_vp, vs .* relerr_vs, geotherms .* relerr_tc1, 
            	crustbase .* relerr_crust, age_uncertainty,
                fill(err_latlong, length(rho)), fill(err_latlong, length(rho)))
            resampled = bsresample(samples, sigma, n, p)
            lats = resampled[:,7]
            longs = resampled[:,8]
        else 
            # Resample! 
            samples = hcat(rho, vp, vs, geotherms, crustbase, ages)
            sigma = hcat(rho .* relerr_rho, vp .* relerr_vp, vs .* relerr_vs, geotherms .* relerr_tc1, 
            	crustbase .* relerr_crust, age_uncertainty)
            resampled = bsresample(samples, sigma, n, p)
        end
        # Save seismic 
        rho = resampled[:,1]
        vp = resampled[:,2]
        vs = resampled[:,3]
        geotherms = resampled[:,4]
        depth = resampled[:,5]
        ages = resampled[:,6]
    end

    if latlong 
        return ((rho, vp, vp ./ vs, geotherms), (depth, ages, lats, longs))
    else
        return ((rho, vp, vp ./ vs, geotherms), (ages))
    end 
end