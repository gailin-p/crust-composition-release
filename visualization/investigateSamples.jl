# For a dataset from commutativityOfAveraging

using ArgParse
using JLD
using MAT
using HDF5
using Statistics
using StatGeochem
using MultivariateStats
using Interpolations
using StatsBase
using Plots; gr();

include("../bin.jl") # TODO locate and use statgeochem version of this

s = ArgParseSettings()
@add_arg_table s begin
    "--data"
        help = "JLD data file saved by commutativityOfAveraging"
        arg_type = String
        required = true
    "--make_figs"
        help = "Figures?"
        arg_type = Bool
        default = false
    "--prefix"
        help = "Prefix for output files, if any (if make_figs or mcign/pca provided)"
        arg_type = String
        default = "../output/average/"
    "--mcign"
        help = "Path to resampled data file, if want to calculate effect on SiO2 estimate"
        arg_type = String
        default = ""
    "--pca"
        help = "Path to JLD with PCA, if want to calculate effect on SiO2 estimate"
        arg_type = String
        default = ""
    "--ign"
        help = "Path to igncn1.mat, if want to calculate effect on SiO2 estimate"
        arg_type = String
        default = "../igncn1.mat"
end
parsed_args = parse_args(ARGS, s)

# order of data: layer, property, run
layers = ["Upper","Middle","Lower"]
props = ["rho","vp","vpvs"] # properties are in this order.

# load data
data = load(parsed_args["data"])
props_of_ave = data["props_of_ave"]
ave_props = data["ave_properties"]
indices = data["indices"]
r = size(indices)[2] # num samples
n = size(indices)[1] # samples per average

# Compare averages
for (l, layer) in enumerate(layers)
    for (p, prop) in enumerate(props)
        ave1 = nanmean(props_of_ave[l,p,:])
        sigma1 = nanstd(props_of_ave[l,p,:])/sqrt(r)
        ave2 = nanmean(ave_props[l,p,:])
        sigma2 = nanstd(ave_props[l,p,:])/sqrt(r)
        Z = (ave1-ave2)/sqrt(sigma1+sigma2)
        println("Z = $Z for $layer, $prop")
    end
end

if parsed_args["make_figs"]
    diffs = props_of_ave - ave_props # layer, property, run (3,3,r)

    for (l, layer) in enumerate(layers)
        for (p, prop) in enumerate(props)
            left = min(minimum(props_of_ave[l,p,:]), minimum(ave_props[l,p,:]))
            right = max(maximum(props_of_ave[l,p,:]), maximum(ave_props[l,p,:]))

            # Stacked histograms: ave of props, props of ave, diffs
            h1 = histogram(props_of_ave[l,p,:], legend=false, xlabel="Property of average", bins=60, xlims=(left,right)) # xlims
            h2 = histogram(ave_props[l,p,:], legend=false, xlabel="Average properties", bins=60, xlims=(left,right))
            h3 = histogram(diffs[l,p,:], legend=false, xlabel="$(prop) difference", bins=60)

            h = plot(h1, h2, h3, layout=(3,1), size=(800,600))

            savefig(h, parsed_args["prefix"]*"compare_hist_$(prop)_$(layer)_n$(n)_r$(r).pdf")
        end
    end
end

### Do inversion to SiO2 as done in runPerplexCrustalStructure
# (TODO move that model to its own module/file so easier to reuse)
# Look at effect of average non-commutativity on SiO2 estimate.
# Compare to actual average SiO2 content.
if (isfile(parsed_args["mcign"]) & isfile(parsed_args["pca"]))
    # What to look at
    iVar = "PC1"
    elem = "SiO2"

    # Import
    println("importing data from $(parsed_args["mcign"]) and $(parsed_args["pca"])")
    ign = matread(parsed_args["ign"])
    models = load(parsed_args["pca"]) # dict layer -> MultivariateStats PCA model
    mcign = Dict{String,Any}()
    for name in ["Age", "Upper_PC1_SiO2", "Middle_PC1_SiO2", "Lower_PC1_SiO2", "SiO2",
        "Calc_Lower_Vp", "Calc_Lower_Rho", "Calc_Lower_VpVs",
        "Calc_Middle_Vp", "Calc_Middle_Rho", "Calc_Middle_VpVs",
        "Calc_Upper_Vp", "Calc_Upper_Rho", "Calc_Upper_VpVs",
        "Calc_Upper_PC1", "Calc_Lower_PC1", "Calc_Middle_PC1"]
        mcign[name] = h5read(parsed_args["mcign"], name)
    end

    # Bin mcign SiO2 by calc PCA to interpolate back from new PCAs to SiO2s
    # (same as in runPerplexCrustalStructure)
    for (l, layer) in enumerate(layers)
        ### SET UP MAPPING

        # Samples for mapping with non-null perplex results
        element_test = map(!isnan, mcign["Calc_" * layer * "_Rho"] .+
                mcign["Calc_" * layer * "_Vp"] .+ mcign["Calc_" * layer * "_VpVs"])

        # Look at independent (seismic) variable of choice
        x = mcign["Calc_" * layer * "_" * iVar]
        xmin = percentile(x, .5)
        xmax = percentile(x,99.5)

        # Bin your element of interest as a function of your seismic variable
        # center (seismic), mean (element of interest eg SiO2), error
        c, m, e = bin(x, mcign[elem][element_test], xmin, xmax,
            length(mcign["SiO2"])/length(ign["SiO2"]), 40)

        # RUN PCA FOR COMMUTATIVITY TEST

        for (name, a) in [("props_of_average", props_of_ave), ("average_properties", ave_props)]
            # Find PCA for each sample. model takes samples like
            #    Array{Float64, 2}(... (3,n)) where rows are rho, vp, vpvs (see runPerplexCrustalStructure)
            Rho_Mean = models[layer * "_Rho_Mean"]
            Rho_Std = models[layer * "_Rho_Std"]
            Vp_Mean = models[layer * "_Vp_Mean"]
            Vp_Std = models[layer * "_Vp_Std"]
            VpVs_Mean = models[layer * "_VpVs_Mean"]
            VpVs_Std = models[layer * "_VpVs_Std"]

            normSampleRho = (a[l,1,:] .- Rho_Mean) ./ Rho_Std;
            normSampleVp = (a[l,2,:] .- Vp_Mean) ./ Vp_Std;
            normSampleVpVs = (a[l,3,:] .- VpVs_Mean) ./ VpVs_Std;

            test = map(!isnan, normSampleRho .+ normSampleVp .+ normSampleVpVs)
            samples = vcat(normSampleRho[test]', normSampleVp[test]', normSampleVpVs[test]')
            sample_pca = transform(models[layer], samples)[1,:]

            # Interpolate (element of interest) values for each (obs seismic) sample in the dataset
            itp = interpolate(m, BSpline(Linear())) # Linear interpolation of bin means
            itp = scale(itp, c) # Scale by bin centers
            itp = extrapolate(itp, Flat()) # outside bin means, just give the upper or lower bin mean. TODO gailin is this the right choice?
            m_itp = itp(sample_pca)

            # Interpolate errors
            itp_e = interpolate(e, BSpline(Linear())) # Linear interpolation of bin means
            itp_e = scale(itp_e, c) # Scale by bin centers
            itp_e = extrapolate(itp_e, Flat()) # TODO gailin - why is it ok to interpolate uncertainty?
            e_itp = itp_e(sample_pca)

            h = histogram(m_itp, xlims=(40,80), xlabel="SiO2 wt%", bins=40)
            savefig(h, parsed_args["prefix"]*"si02-comm-n$n-r$r-$(layer)-$(name).pdf")
            println("n=$n average wt% SiO2 for $(layer) $(name) = $(nanmean(m_itp)), with uncertainty $(nanmean(e_itp)))")
        end
    end

    # Compare to actual SiO2 mean over entire dataset
    println("Mean SiO2 over all samples is $(mean(ign["SiO2"]))")

    # TODO compare to mean SiO2 for each set of samples
else
    println("Not running SiO2 analysis")
end
