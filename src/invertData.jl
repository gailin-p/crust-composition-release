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
  FRSeismic()
Same idea as getAllSeismic, but return fountain/rudnick data. Only returns vp.
layer: 6, 7, or 8 (upper, middle, lower)
"""
function FRSeismic(layer::Integer; model::Integer=2, n::Integer=400, resample::Bool=true)
  dat, header = isfile("Fountain-Rudnick.csv") ? readdlm("Fountain-Rudnick.csv", ',', header=true) :
    readdlm("../Fountain-Rudnick.csv", ',', header=true)
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
  getPerplexTestSeismic()

Fake earth builder, using Perple_X modeling to model earth

"""
function getPerplexTestSeismic(n, layer; nsamples=100, target=[66.6, 63.5, 53.4], bonus_info=false)
  # Load pre-resampled ign data
  layer = layer - 5
  f = isfile("resources/bsr_ignmajors_1.csv") ? "resources/bsr_ignmajors_1.csv" : "../resources/bsr_ignmajors_1.csv"
  dat, h = readdlm(f, ',', header=true)
  sio2 = dat[:,2]
  indices = Int.(dat[:,1])


  # Iterative method to find set of nsamples rocks that match target composition at this layer
  s = sample(indices,nsamples, replace=false)
  snew = zeros(Int, length(s))
  snew .= s
  best = abs(mean(sio2[s]) - target[layer])
  while best > 1/20
      to_replace = sample(1:nsamples)
      new = sample(indices)
      snew[to_replace] = new
      if abs(mean(sio2[snew]) - target[layer]) < best
          s .= snew
          best = abs(mean(sio2[snew]) - target[layer])
          #println("new best $(mean(sio2[snew]))")
      end
      snew .= s # reset
  end
  s = Int.(s)
  # Major element comp
  comps = dat[s, 2:11]

  ## Get crust properties
  crust = crustDistribution.getCrustParams(n, uncertain=true) # depth to 550, base of upper, middle, lower
  geotherm = crust[:,1]
  boundaries = hcat(fill(0, n), crust[:,2:end]) # first boundary is 0
  depths = ((boundaries[:,layer+1] .- boundaries[:,layer]) ./ 2) .+ boundaries[:,layer] # middle of layer

  ### Assign composition to each sample location
  assignment = sample(1:nsamples, n)
  ### Log assigned compositions for debugging
  @info "Assigned comps" comps[assignment,:]

  # run Perple_X for each sample, then adjust to props at each location assigned that sample
  st = SeismicTransform()
  rho = zeros(n); vp = zeros(n); vs = zeros(n);

  failed = zeros(0)
  succeeded = zeros(0)
  @showprogress "Building fake earth" for sample in 1:nsamples
    ####### Shared
    dtdz = 550.0/nanmean(geotherm) # sampling geotherm
    for t in 1:10 # try up to 10 formation depths and pressures
      formation_depth, formation_dtdz = crustDistribution.getFormationParams(nanmean(depths), dtdz)

      P_range = [formation_depth*(9/10)*dpdz, formation_depth*(11/10)*dpdz] # we only need a small range around formation t, p
      # Run perplex
      perplex_configure_geotherm(DEFAULT_PERPLEX, DEFAULT_SCRATCH, comps[sample,:], PERPLEX_COMPOSITION_ELTS,
                    P_range, 273.15, formation_dtdz/dpdz, dataset=DEFAULT_DATASET, solution_phases=SOLUTIONS,
                    excludes=FLUID_ENDMEMBERS, index=1, npoints=NPOINTS)
      point = perplex_query_point(DEFAULT_PERPLEX, DEFAULT_SCRATCH, formation_depth*dpdz, index=1)
      try
        properties = get_system_props(point)
        endmembers = parse_perplex_point(point)

        for (j, loc) in enumerate(assignment)
          if assignment[loc] == sample
            P = dpdz*depths[j] + 280 # add surface pressure (bar)
            T = geotherm[j]*depths[j] + 273.15 # add surface temp (K)
            rho[j], vp[j], vs[j] = get_seismic(T, P, properties, endmembers, st)
          end
        end
        push!(succeeded,comps[sample,1])
        break # successful
      catch e
        #if (isa(e, ParsePerplexError) | isa(e, SeismicError))
        println("\r\n\r\nCannot process sample due to \r\n $e")
        println("failed sample had $(comps[sample,1]) % sio2, try $t")
        if t == 10
          push!(failed,comps[sample,1])
        end
        continue
        #else
        # throw(e)
        #end
      end
    end
  end

  println("layer $layer has $(length(failed)) failed samples, mean sio2 $(mean(failed))")

  rho[isnan.(rho)] .= nanmean(rho)
  vp[isnan.(vp)] .= nanmean(vp)
  vs[isnan.(vs)] .= nanmean(vs)

  # Cracking in upper crust
  if layer == 1
	println("Applying cracking to upper crust of simulated earth")
    liquid_weights = [.5, .5, 0, 0] # fraction dry, water, magma, no cracking
    crack_porosity = .007 * .05 # total porosity * fraction crack
    pore_porosity = .007 * (1-.05) # total porosity * fraction pore

    profiles = Array{Crack,1}([random_cracking(crack_porosity, pore_porosity, liquid_weights) for i in 1:n])
    vpvs = vp ./ vs
    rho, vp, vpvs = apply_cracking!(rho, vp, vpvs, profiles)
    return rho, vp, vpvs, geotherm
  end

  if bonus_info
    return (rho, vp, vp ./ vs, geotherm), (s, failed, succeeded)
  else
    return (rho, vp, vp ./ vs, geotherm)
  end
end


# Return (rho, vp, vp/vs, geotherm)
# Uncertainty is fraction of Crust1.0 std for this layer to apply as std of gaussian noise.
function getTestSeismic(n, layer, uncertainty, samples=["D95-13", "D95-16", "HL3"])
    if !(layer in [6,7,8])
        throw(ArgumentError("Layer must be 6, 7, or 8 (crysteline Crust1 layers)"))
    end

    if isfile("resources/dabie/kern_dabie_comp.csv")
      compdat, comph = readdlm("resources/dabie/kern_dabie_comp.csv", ',', header=true)
      seismicdat, seismich = readdlm("resources/dabie/kern_dabie_seismic.csv", ',', header=true);
    elseif isfile("../resources/dabie/kern_dabie_comp.csv")
      compdat, comph = readdlm("../resources/dabie/kern_dabie_comp.csv", ',', header=true)
      seismicdat, seismich = readdlm("../resources/dabie/kern_dabie_seismic.csv", ',', header=true);
    else
      throw(Error("Cannot find file kern_dabie_comp.csv in data dir"))
    end
    comph = comph[:]
    seismich = seismich[:]

    # Gather data...
    crust = crustDistribution.getCrustParams(n, uncertain=true) # depth to 550, base of upper, middle, lower
    geotherm = crust[:,1]
    boundaries = hcat(fill(0, n), crust[:,2:end]) # first boundary is 0
    layer = layer - 5 # now layer of interest is 1, 2, or 3
    depths = ((boundaries[:,layer+1] .- boundaries[:,layer]) ./ 2) .+ boundaries[:,layer] # middle of layer

    # Sample params
    if compdat[:,1] != seismicdat[:,1]
        error("Sample names do not line up for Dabie composition and seismic files.")
    end
    samplei = findfirst(isequal(samples[layer]), compdat[:,1])
    rhoj = findfirst(isequal("rho (g/cm^3)"), seismich)

    vp0 = seismicdat[samplei, findfirst(isequal("Vp0"), seismich)]
    vs0 = seismicdat[samplei, findfirst(isequal("Vs0"), seismich)]
    dvpdp = seismicdat[samplei, findfirst(isequal("dVp=dP"), seismich)]*(10^-4)/10 # MPa to bar
    dvsdp = seismicdat[samplei, findfirst(isequal("dVs=dP"), seismich)]*(10^-4)/10 # MPa to bar
    dvpdt = seismicdat[samplei, findfirst(isequal("dVp=dT"), seismich)]*(-10^-4)
    dvsdt = seismicdat[samplei, findfirst(isequal("dVs=dT"), seismich)]*(-10^-4)

    # For each crust point, estimate vp and vs at depths
    rho = fill(seismicdat[samplei, rhoj]*1000, n)

    dtdz = 550 ./ geotherm

    vp = vp0 .+ (dvpdp .* (dpdz .* depths)) .+ (dvpdt .* (dtdz .* depths))
    vs = vs0 .+ (dvsdp .* (dpdz .* depths)) .+ (dvsdt .* (dtdz .* depths))

    # Use uncertainty as portion of crust1 std for comparability
    (vp_c, vs_c, rho_c) = find_crust1_seismic(crustDistribution.all_lats, crustDistribution.all_longs, layer+5)
    rho = rho .+ uncertainty*std(rho_c)*randn(n)
    vp = vp .+ uncertainty*std(vp_c)*randn(n)
    vs = vs .+ uncertainty*std(vs_c)*randn(n)

    return (rho, vp, vp ./ vs, geotherm)

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

    elseif dataSrc == "Shen"
      # Get data
      # Dims: longitude latitude depth. NaN-value = 9999.0
      vp_us = Array{Float64}(ncread("resources/shen.US.2016.nc", "vp")) .|> x -> x==9999.0 ? NaN : x
      vs_us = Array{Float64}(ncread("resources/shen.US.2016.nc", "vsv")) .|> x -> x==9999.0 ? NaN : x
      rho_us = Array{Float64}(ncread("resources/shen.US.2016.nc", "rho")) .|> x -> x==9999.0 ? NaN : x
      rho_us = rho_us .* 1000 # unit conversion

      # values corresponding to each axis of the above data
      depth_us = Array{Float64}(ncread("resources/shen.US.2016.nc", "depth"))
      latitude_us = Array{Float64}(ncread("resources/shen.US.2016.nc", "latitude"))
      # values range, bizarely, from 235 to 295. This puts it in the right range but idk if correct.
      longitude_us = Array{Float64}(ncread("resources/shen.US.2016.nc", "longitude") .- 360)

      # Get data as 1d array
      lats = repeat(latitude_us, inner=length(longitude_us)) # repeat each lat for length of longs
      longs = repeat(longitude_us, outer=length(latitude_us)) # repeat sequence of lats for each lat
      vp = fill(NaN, length(lats))
      vs = fill(NaN, length(lats))
      rho = fill(NaN, length(lats))

      layer_tops = -1 .* find_crust1_base(lats, longs, layer-1) # for each lat/long, top of layer
      layer_bottoms = -1 .* find_crust1_base(lats, longs, layer) # for each lat/long, bottom of layer
      i = 1
      for lat in 1:length(latitude_us)
        for long in 1:length(longitude_us)
          # Depth of middle of this layer at this lat/long
          depth = layer_tops[i] + (layer_bottoms[i] - layer_tops[i])/2
          #println("depth of center of layer $(depth)")

          # Corresponding index of depth availible from Shen seismic data
          shen_depth = findfirst(x->x>=depth, depth_us)
          if rand() >= .5
            shen_depth -= 1 # half the time, use the shallower option
          end

          vp[i] = vp_us[long,lat,shen_depth]
          vs[i] = vs_us[long,lat,shen_depth]
          rho[i] = rho_us[long,lat,shen_depth]

          i += 1
        end
      end

      geotherms = find_tc1_crust(lats, longs) # Geotherm
      crustbase = -1 .* find_crust1_base(lats, longs, 8) # base of lower crust

      # TODO nicer to check before doing calculations on all lats longs but I'm lazy
      cont = map(c -> continents[c], find_geolcont(lats,longs))
      good = .~ (isnan.(vp) .| isnan.(rho) .| isnan.(vs) .| isnan.(geotherms) .| (cont .== "NA"))
      println("$(sum(good)) good values from Shen et al.")
      vp = vp[good]
      vs = vs[good]
      rho = rho[good]
      lats = lats[good]
      longs = longs[good]
      geotherms = geotherms[good]
      crustbase=crustbase[good]

    elseif (dataSrc == "Dabie") | (dataSrc == "DabieRG") | (dataSrc == "TestRG") ## Special case: "fake earth" where we know target comp.
      if dataSrc == "Dabie"
        rho, vp, vpvs, geotherms = getTestSeismic(n, layer, dataUncertainty)
      elseif dataSrc == "DabieRG"
        rho, vp, vpvs, geotherms = getTestSeismic(n, layer, dataUncertainty, ["D95-44", "HT4", "D95-11"])
      else  # testRG
        rho, vp, vpvs, geotherms = getPerplexTestSeismic(n, layer)
      end
      lats = fill(NaN, n)
      longs = fill(NaN, n)
      ages = fill(NaN, n)
      crustbase = fill(NaN, n)
      if latlong
        return ((rho, vp, vpvs, geotherms), (crustbase, ages, lats, longs))
      else
          return ((rho, vp, vpvs, geotherms), (ages))
      end

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
