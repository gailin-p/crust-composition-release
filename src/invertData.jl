"""
Data from Crust1.0 to invert.

Allows choice of age models.

Requires user to load config.jl and crustDistribution.jl
"""


using MAT
using NetCDF
using StatGeochem
using StatsBase
using Statistics
using Test
using Random
using ProgressMeter: @showprogress

# one std relative error. If running from jupyter notebook, need relative path
uncertainty_dat = matread(IGN_FILE)["err2srel"]

# Now using std as error, see getAllSeismic
# relerr_rho = uncertainty_dat["Rho"]/2
# relerr_vp = uncertainty_dat["Vp"]/2
# relerr_vs = uncertainty_dat["Vs"]/2
# relerr_rho = .05 # Huang et al error estimate: 1 std = 5% err
# relerr_vp = .05
# relerr_vs = .05
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
  getTestSeismic()
Build a test earth directly from model seismic properties: so model contains superset
of rocks in the test earth
"""
function getTestSeismic(n, models::ModelCollection, target::Float64, layer::String; nsamples=595)
    # compositions to choose from
    comps = models.models[layer][1].comp
    indices = Array(1:size(comps,1))
    sio2 = comps[:,2]

    # choose comp samples from ign of model 1 that average to target
	# while ensuring all have non-nan seismic props in all geotherm bins
	s = fill(0, nsamples); i = 1
	while i <= nsamples
		new = sample(indices)
		# Are seismic props in all geotherm bins non-nan for this sample?
		seismic_ok = true
		for m in models.models[layer]
			if isnan(sum(m.seismic[new,:]))
			  seismic_ok = false
			  #println("discard sample $new with props $(m.seismic[new,:])")
			end
		end
		if seismic_ok
			s[i] = new
			#println("assignmed $i")
			i += 1
		end
	end

    snew = zeros(Int, length(s))
    snew .= s
    best = abs(mean(sio2[s]) - target)
    while best > 1/100
      to_replace = sample(1:nsamples)
      new = sample(indices)
      snew[to_replace] = new

	  # Are seismic props in all geotherm bins non-nan for this sample?
	  seismic_ok = true
	  for m in models.models[layer]
		  if isnan(sum(m.seismic[new,:]))
			  seismic_ok = false
		  end
	  end

	  # Does the error decrease with inclusion of this sample?
      if (abs(mean(sio2[snew]) - target) < best) & seismic_ok
          s .= snew
          best = abs(mean(sio2[snew]) - target)
          #println("new best $(mean(sio2[snew]))")
      end
      snew .= s # reset
    end
    s = Int.(s)

    # n geotherm samples
    crust = crustDistribution.getCrustParams(n, uncertain=true) # depth to 550, base of upper, middle, lower
    geotherm = crust[:,1]
    bins = [searchsortedlast(Array(models.bins), g) for g in geotherm]
	bins[bins .== 11] .= 10 # Because geotherms are noisily assigned, some might be smaller or larger than actual geotherm bins. Bound them.
	bins[bins .== 0] .= 1

    # assign comp to each of n
    assignment = sample(s, n)
	@info "Assigned comps" comps[assignment,2:11]

    # for each of n, choose the seismic props from the model with the correct geotherm bin
    seismic = zeros(n,3)
    for i in 1:n
        seismic[i,:] .= models.models[layer][bins[i]].seismic[assignment[i],2:end]
    end
    return (seismic[:,1], seismic[:,2], seismic[:,3], geotherm), (fill(NaN,n), fill(NaN,n),fill(NaN,n),fill(NaN,n))
end


"""
    getAllSeismic()
Return a touple of weights and lists of seismic values (rho, vp, vp/vs, geotherm)
From Crust1.0
Optionally resample (using weights corresponding to area of each 1x1 degree square).
Latlong info optionally returned.
By default, resample with Tc1 ages
Optionally add a systematic bias to Crust1.0 seismic data after resampling.

TODO THIS IS TERRIBLE SPAGHETTI CODE, REFACTOR ME !
"""
function getAllSeismic(layer::Integer;
  ageModel::AgeModel=Tc1Age(), n::Integer=50000, resample::Bool=true, systematic::Bool=false,
  latlong::Bool=false, dataSrc::String="Crust1.0", dataUncertainty::Float64=1.0)
    if !(layer in [6,7,8])
        throw(ArgumentError("Layer must be 6, 7, or 8 (crysteline Crust1 layers)"))
    end

   if dataSrc == "Crust1.0"
      # Gather data...
      lats = crustDistribution.all_lats # Lats where both tc1 and crust1.0 defined
      longs = crustDistribution.all_longs
      geotherms = crustDistribution.depth[:,1]
      crustbase = crustDistribution.depth[:,4]

      (vp, vs, rho) = find_crust1_seismic(lats, longs, layer)
  elseif dataSrc == "Spiral"
	  return getSpiralSeismic(layer, ageModel, n, latlong)

    else
      throw("Unrecognized data source")
    end

    ages, age_uncertainty = sampleAge(lats, longs, ageModel)

    err_rho = std(rho) * dataUncertainty
    err_vp = std(vp) * dataUncertainty
    err_vs = std(vs) * dataUncertainty

    if resample
        k = crustDistribution.latLongWeight.(lats)
        # Probability of keeping a given data point when sampling:
        # We want to select roughly one-fith of the full dataset in each re-sample,
        # which means an average resampling probability <p> of about 0.2
        #p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
        p = k .* .2/mean(k)

        if latlong
            samples = hcat(rho, vp, vs, geotherms, crustbase, ages, lats, longs)
            sigma = hcat(fill(err_rho, length(rho)), fill(err_vp, length(vp)), fill(err_vs, length(vs)), geotherms .* relerr_tc1,
            	crustbase .* relerr_crust, age_uncertainty,
                fill(err_latlong, length(rho)), fill(err_latlong, length(rho)))
            resampled = bsresample(samples, sigma, n, p)
            lats = resampled[:,7]
            longs = resampled[:,8]
        else
            # Resample!
            samples = hcat(rho, vp, vs, geotherms, crustbase, ages)
            sigma = hcat(fill(err_rho, length(rho)), fill(err_vp, length(vp)), fill(err_vs, length(vs)), geotherms .* relerr_tc1,
            	crustbase .* relerr_crust, age_uncertainty)
            resampled = bsresample(samples, sigma, n, p)
        end
        # Save seismic
        rho = resampled[:,1]
        vp = resampled[:,2]
        vs = resampled[:,3]
        geotherms = resampled[:,4]
        crustbase = resampled[:,5]
        ages = resampled[:,6]
    end

    if systematic
      rho = rho .+ err_rho*randn()
      vp = vp .+ err_vp*randn()
      vs = vs .+ err_vs*randn()
    end

    if latlong
        return ((rho, vp, vp ./ vs, geotherms), (crustbase, ages, lats, longs))
    else
        return ((rho, vp, vp ./ vs, geotherms), (ages))
    end
end


function getSpiralSeismic(layer::Integer, ageModel::AgeModel=Tc1Age(),
	n::Integer=50000, latlong::Bool=false)

	## Sample n locations on continents
	idx = sample(Array(1:length(crustDistribution.all_lats)),
		Weights(crustDistribution.latLongWeight.(crustDistribution.all_lats)), n)

	lats = crustDistribution.all_lats[idx]
	longs = crustDistribution.all_longs[idx]
	geotherms = crustDistribution.depth[idx,1]

	# Load files into maps for easy lat/long lookup
	if layer == 6
		file = "resources/spiral/SPiRaL_1.4.Interpolated.Surface.12.Crust_12_Bottom_of_upper_crust.txt"
		lith_file = "resources/litho1.0/lith_layer6.csv"
	elseif layer == 7
		file = "resources/spiral/SPiRaL_1.4.Interpolated.Surface.14.Crust_14_Bottom_of_middle_crust.txt"
		lith_file = "resources/litho1.0/lith_layer7.csv"
	elseif layer == 8
		file = "resources/spiral/SPiRaL_1.4.Interpolated.Surface.16.Crust_16_Bottom_of_lower_crust.txt"
		lith_file = "resources/litho1.0/lith_layer8.csv"
	else
		throw("Unrecognized layer option")
	end
	dat = readdlm(file)
	indices = readdlm("resources/spiral/SPiRaL_1.4.Interpolated.Coordinates.txt")

	lookup = Dict{Tuple{Float64, Float64},Tuple{Float64, Float64, Float64}}()
	for i in 1:size(dat,1)
		# mean vp, mean vs, base of layer depth
		lookup[(indices[i,1],indices[i,2])] = ((dat[i,3]+dat[i,4])/2, (dat[i,5]+dat[i,6])/2, dat[2])
	end

	vp = [lookup[(lats[i], longs[i])][1] for i in 1:length(lats)]
	vs = [lookup[(lats[i], longs[i])][2] for i in 1:length(lats)]
	crustbase = [lookup[(lats[i], longs[i])][3] for i in 1:length(lats)]

	# rho from crust1
	rho_dat, _ = readdlm(lith_file, ',', header=true)
	rho_lookup = Dict{Tuple{Float64, Float64},Float64}() # construct lookup dict
	for i in 1:size(rho_dat,1)
		rho_lookup[(rho_dat[i,1],rho_dat[i,2])] = rho_dat[i,5]
	end
	rho = [rho_lookup[(lats[i], longs[i])] for i in 1:length(lats)]

	ages, age_uncertainty = sampleAge(lats, longs, ageModel)

	if latlong
        return ((rho, vp, vp ./ vs, geotherms), (crustbase, ages, lats, longs))
    else
        return ((rho, vp, vp ./ vs, geotherms), (ages))
    end
end
