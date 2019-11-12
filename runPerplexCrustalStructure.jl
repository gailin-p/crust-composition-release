## --- Load requisite Julia packages (may need to install first)

    # Load (and install if neccesary) the StatGeochem package which has the resampling functions we"ll want
    try
        using StatGeochem
    catch
        using Pkg
        Pkg.clone("https://github.com/brenhinkeller/StatGeochem.jl")
        using StatGeochem
    end

    using Statistics, DelimitedFiles, MAT, StatsBase, Interpolations
    using ProgressMeter: @showprogress
    using Plots; gr();
    using Logging
    using MultivariateStats
    using Random
    using HDF5

    include("bin.jl")

# ## --- Export out element compositions to run runPerplexBatchVp on cluster
#
#     # Read igneous whole-rock geochemistry from .mat file
#     ign = matread("igncn1.mat")
#     ign["elements"] = string.(ign["elements"]) # Since .mat apparently can"t tell the difference between "B" and "b"
#
#     # List of elements we want to export (primarily out element oxides)
#     elements = ["index","SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2","tc1Crust"]
#
#     # Copy over just the elements we want to export
#     out = Dict()
#     for e in elements
#         out[e] = ign[e]
#     end
#
#     # Start with all Fe as FeO
#     out["FeO"] = feoconversion.(ign["FeO"], ign["Fe2O3"], ign["FeOT"], ign["Fe2O3T"])
#
#     # Set undefined H2O to the average
#     out["H2O_Total"][isnan.(out["H2O_Total"])] .= nanmean(out["H2O_Total"])
#
#     # Set undefined CO2 to the average
#     out["CO2"][isnan.(out["CO2"])] .= nanmean(out["CO2"])
#
#     # Create 2d array of out element data to export
#     outtable = Array{Float64,2}(undef, length(out["SiO2"]), length(elements))
#     for i = 1:length(elements)
#         outtable[:,i] = out[elements[i]]
#     end
#
#     # Reject samples with missing data
#     t = .~ any(isnan.(outtable), dims=2)
#
#     # Reject samples with suspicious anhydrous normalizations
#     anhydrousnorm = sum(outtable[:,2:9], dims=2)
#     t = t .& (anhydrousnorm .< 101) .& (anhydrousnorm .> 90)
#
#     # Write accepted samples to file
#     if isfile("ignmajors.csv")
#         @warn "ignmajors.csv already exists, was not overwritten!"
#     else
#         writedlm("ignmajors.csv", round.(outtable[t[:],:], digits=5), ",")
#     end
#
# # --- Run Perplex on cluster for each sample
#
#     Use runPerplexBatchVp and runProduction.pbs
#     Get "PerplexResults.log" output file.
#
# --- Import results from PerplexResults.log

    # Read igneous whole-rock geochemistry from .mat file
    ign = matread("igncn1.mat")
    ign["elements"] = string.(ign["elements"]) # Since .mat apparently can"t tell the difference between "B" and "b"

    # Read PerplexResults.log
    perplexresults = importdataset("data/runTestNew.small.log", '\t')

    # Index must be an integer
    perplexresults["index"] = Int.(perplexresults["index"])

    # Initialize calculated variables with NaNs
    for e in ["Calc_Upper_Rho", "Calc_Upper_Vp", "Calc_Upper_VpVs", "Calc_Middle_Rho", "Calc_Middle_Vp", "Calc_Middle_VpVs", "Calc_Lower_Rho", "Calc_Lower_Vp", "Calc_Lower_VpVs",]
        ign[e] = fill(NaN, length(ign["SiO2"]))
        ign["err2srel"][e] = fill(NaN, length(ign["SiO2"]))
    end

    # Pressure gradient (bar / km) -- we"ll need this soon
    dpdz = 2818 * 9.8 / 10^5 *10^3

    # Assign calculated values
    @showprogress 1 "Importing PerpleX results: " for i in unique(perplexresults["index"])
        # Find calculated seismic properties for sample with index i
        # GAILIN - added filter here b/c these are very high resolution (1193 points per eartchem sample)
        all = 1:1:length(perplexresults["index"])
        t = (perplexresults["index"] .== i) .& ((all .% 10) .== 0)

        # Isolate the calculated seismic properties for sample i
        rho = perplexresults["rho"][t]
        vp =  perplexresults["Vp"][t]
        vpvs =  perplexresults["Vp/Vs"][t]
        pressure = perplexresults["P(bar)"][t]

        # Find the row containing the sample with index i
        row = findfirst(ign["index"] .== i)

        # Upper crust pressure range
        upper = pressure .<= ign["Upper_Crust"][row] * dpdz
        ign["Calc_Upper_Rho"][row] = nanmean(rho[upper])
        ign["Calc_Upper_Vp"][row] = nanmean(vp[upper])
        ign["Calc_Upper_VpVs"][row] = nanmean(vpvs[upper])
        ign["err2srel"]["Calc_Upper_Rho"][row] = nanstd(rho[upper]) ./ ign["Calc_Upper_Rho"][row] * 2
        ign["err2srel"]["Calc_Upper_Vp"][row] = nanstd(vp[upper]) ./ ign["Calc_Upper_Vp"][row] * 2
        ign["err2srel"]["Calc_Upper_VpVs"][row] = nanstd(vpvs[upper]) ./ ign["Calc_Upper_VpVs"][row] * 2


        # Middle crust pressure range
        middle = (pressure .> ign["Upper_Crust"][row] * dpdz) .& (pressure .<= (ign["Upper_Crust"][row] + ign["Middle_Crust"][row]) * dpdz)
        ign["Calc_Middle_Rho"][row] = nanmean(rho[middle])
        ign["Calc_Middle_Vp"][row] = nanmean(vp[middle])
        ign["Calc_Middle_VpVs"][row] = nanmean(vpvs[middle])
        ign["err2srel"]["Calc_Middle_Rho"][row] = nanstd(rho[middle]) ./ ign["Calc_Middle_Rho"][row] * 2
        ign["err2srel"]["Calc_Middle_Vp"][row] = nanstd(vp[middle]) ./ ign["Calc_Middle_Vp"][row] * 2
        ign["err2srel"]["Calc_Middle_VpVs"][row] = nanstd(vpvs[middle]) ./ ign["Calc_Middle_VpVs"][row] * 2

        # Lower crust pressure range
        lower = (pressure .> (ign["Upper_Crust"][row] + ign["Middle_Crust"][row]) * dpdz) .& (pressure .<= ign["Crust"][row] * dpdz)
        ign["Calc_Lower_Rho"][row] = nanmean(rho[lower])
        ign["Calc_Lower_Vp"][row] = nanmean(vp[lower])
        ign["Calc_Lower_VpVs"][row] = nanmean(vpvs[lower])
        ign["err2srel"]["Calc_Lower_Rho"][row] = nanstd(rho[lower]) ./ ign["Calc_Lower_Rho"][row] * 2
        ign["err2srel"]["Calc_Lower_Vp"][row] = nanstd(vp[lower]) ./ ign["Calc_Lower_Vp"][row] * 2
        ign["err2srel"]["Calc_Lower_VpVs"][row] = nanstd(vpvs[lower]) ./ ign["Calc_Lower_VpVs"][row] * 2
    end

## --- Bootstrap resampling

    # Calculate spatiotemporal sample density factor ("inverse weight")
    k = invweight(ign["Latitude"], ign["Longitude"], ign["Age"])

    # Probability of keeping a given data point when sampling:
    # We want to select roughly one-fith of the full dataset in each re-sample,
    # which means an average resampling probability <p> of about 0.2
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)

    # Absolute uncertainties for each variable:
    # Start with default 2-sigma relative uncertainties from err2srel
    for e in ign["elements"]
        ign["err"][e] = ign[e] .* (ign["err2srel"][e] / 2)
    end

    # Special cases: Latitude & Longitude
    ign["err"]["Latitude"] = ign["Loc_Prec"]
    ign["err"]["Longitude"] = ign["Loc_Prec"]

    # Special cases: Age
    ign["err"]["Age"] = (ign["Age_Max"]-ign["Age_Min"])/2;
    # Find points with < 50 Ma absolute uncertainty
    t = (ign["err"]["Age"] .< 50) .| isnan.(ign["err"]["Age"])
    # Set 50 Ma minimum age uncertainty (1-sigma)
    ign["err"]["Age"][t] .= 50;

## --  Resample a few hundred times (all elements at once)

    elements = ["index", "Latitude", "Longitude", "Elevation", "Age", "SiO2", "TiO2", "Al2O3", "Fe2O3", "Fe2O3T", "FeO", "FeOT", "MgO", "CaO", "Na2O", "K2O", "P2O5", "MnO", "H2O_Total", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Li", "Be", "B", "C", "CO2", "F", "Cl", "Sc", "Ti", "V", "Cr", "Co", "Ni", "Cu", "Zn", "Ga", "Zr", "Os", "Rb", "Bi", "Hg", "Ba", "Y", "Pb", "Te", "Nb", "Sr87_Sr86", "Tl", "Pt", "Sn", "Cd", "As", "Pd", "Sr", "Se", "S", "Au", "Ta", "Mo", "U", "Cs", "Sb", "Ag", "W", "Th", "Re", "Hf", "Ir", "tc1Lith", "tc1Crust", "Crust", "Vp", "Vs", "Rho", "Upper_Crust", "Upper_Vp", "Upper_Vs", "Upper_Rho", "Middle_Crust", "Middle_Vp", "Middle_Vs", "Middle_Rho", "Lower_Crust", "Lower_Vp", "Lower_Vs", "Lower_Rho", "Freeair", "Bouger", "Eustar", "Calc_Upper_Vp", "Calc_Upper_VpVs", "Calc_Upper_Rho", "Calc_Middle_Vp", "Calc_Middle_VpVs", "Calc_Middle_Rho", "Calc_Lower_Vp", "Calc_Lower_VpVs", "Calc_Lower_Rho"]
    nrows = 10^7 # Resample up to...
    mcign = bsresample(ign, nrows, elements, p)

## --- Some data quality checks

@info "$(length(findall(isnan, mcign["Calc_Upper_Vp"]))) NaN values out of
    $(length(mcign["Calc_Upper_Vp"])) total for calculated, resampled property Calc_Upper_Vp"

@info "$(length(findall(isnan, mcign["Upper_Vp"]))) NaN values out of
    $(length(mcign["Upper_Vp"])) resampled property Upper_Vp"


## --- PCA
# GAILIN - PCA of the seismic properties given by perple_x

features = ["Rho", "Vp", "VpVs"]
layers = ["Upper", "Middle", "Lower"]

@showprogress "Running PCA for upper, middle, lower crust: " for crust in layers
#crust = "Upper"

    Rho_Mean = nanmean(mcign["Calc_" * crust * "_Rho"]);
    Rho_Std = nanstd(mcign["Calc_" * crust * "_Rho"]);
    Vp_Mean = nanmean(mcign["Calc_" * crust * "_Vp"]);
    Vp_Std = nanstd(mcign["Calc_" * crust * "_Vp"]);
    VpVs_Mean = nanmean(mcign["Calc_" * crust * "_VpVs"]);
    VpVs_Std = nanstd(mcign["Calc_" * crust * "_VpVs"]);

    CalcRhoNorm = (mcign["Calc_" * crust * "_Rho"] .- Rho_Mean) ./ Rho_Std;
    CalcVpNorm = (mcign["Calc_" * crust * "_Vp"] .- Vp_Mean) ./ Vp_Std;
    CalcVpVsNorm = (mcign["Calc_" * crust * "_VpVs"] .- VpVs_Mean) ./ VpVs_Std;

    test = map(!isnan, CalcRhoNorm .+ CalcVpNorm .+ CalcVpVsNorm)

    samples = vcat(CalcRhoNorm[test]', CalcVpNorm[test]', CalcVpVsNorm[test]')

    # PCA
    model = fit(PCA, samples)
    mcign["Calc_" * crust * "_PC1"] = transform(model, samples)[1,:]

    # Add VpVs, if necessary
    if ~(crust*"_VpVs" in keys(mcign))
        mcign[crust *"_VpVs"] = mcign[crust * "_Vp"] ./ mcign[crust * "_Vs"]
    end

    # Make sure density is in kg/m3
    if nanmean(mcign[crust * "_Rho"]) < 1000
        mcign[crust * "_Rho"] = mcign[crust * "_Rho"] .* 1000
    end

    ObsRhoNorm = (mcign[crust * "_Rho"] .- Rho_Mean) ./ Rho_Std
    ObsVpNorm = (mcign[crust * "_Vp"] .- Vp_Mean) ./ Vp_Std
    ObsVpVsNorm = (mcign[crust * "_VpVs"] .- VpVs_Mean) ./ VpVs_Std

    obs = vcat(ObsRhoNorm', ObsVpNorm', ObsVpVsNorm')

    # Do PCA for samples resampled from crust1
    mcign[crust * "_PC1"] = transform(model, obs)[1,:] # Save PC1

end

## --- Estimate compositions
nbins = 8
tmin = 0
tmax = 4000
iVar = "PC1"
elem = "SiO2"

# For each crust type
#@showprogress "Running composition: " for crust in layers
p = plot();
for crust in layers
    mcign[crust * '_' * iVar * '_' * elem] = fill(NaN, length(mcign["Age"]));

    # For each age bin
    rt = tmin:(tmax-tmin)/nbins:tmax
    for j = 1:nbins
        test = (mcign["Age"] .> rt[j]) .& (mcign["Age"] .< rt[j+1])
        element_test = map(!isnan, mcign["Calc_" * crust * "_Rho"] .+
                mcign["Calc_" * crust * "_Vp"] .+ mcign["Calc_" * crust * "_VpVs"])

        # Look at your independent (seismic) variable of choice
        x = mcign["Calc_" * crust * "_" * iVar]
        xmin = percentile(x, .5)
        xmax = percentile(x,99.5)

        # Bin your element of interest as a function of your seismic variable
        # center (seismic), mean (element of interest eg SiO2), error
        c, m, e = bin(x, mcign[elem][element_test], xmin, xmax,
            length(mcign["SiO2"])/length(ign["SiO2"]), 40)

        # Interpolate (element of interest) values for each (obs seismic) sample in the dataset
        itp = interpolate(m, BSpline(Linear())) # Linear interpolation of bin means
        sitp = scale(itp, c) # Scale by bin centers
        esitp = extrapolate(sitp, Flat()) # outside bin means, just give the upper or lower bin mean. TODO gailin is this the right choice?
        m_itp = esitp(mcign[crust * "_" * iVar][test])

        # Interpolate errors
        itp_e = interpolate(e, BSpline(Linear())) # Linear interpolation of bin means
        sitp_e = scale(itp_e, c) # Scale by bin centers
        esitp_e = extrapolate(sitp_e, Flat()) # TODO gailin - why is it ok to interpolate uncertainty?
        e_itp = esitp_e(mcign[crust * "_" * iVar][test])

        # Propagate uncertainty
        # TODO gailin I don't really understand this
        mcign[crust * "_" * iVar * "_" * elem][test] = m_itp .+ (rand!(fill(NaN, length(m_itp))) .* e_itp) # random errors scaled by interpolated errors?
    end

    #% Plot the temporal trend for this crust type
    age_centers, elt_means, elt_errors = bin(mcign["Age"], mcign[crust * "_" * iVar * "_" * elem],
            tmin, tmax, length(mcign["SiO2"])/length(ign["SiO2"]), nbins)
    plot!(age_centers, elt_means, yerror=elt_errors, label=crust);
end

# Plot exposed
age_centers, elt_means, elt_errors = bin(mcign["Age"], mcign[elem],
        tmin, tmax, length(mcign["SiO2"])/length(ign["SiO2"]), nbins)
plot!(age_centers, elt_means, yerror=elt_errors, label="Exposed");

plot!(xlabel="Age");
plot!(ylabel="SiO2 composition");
plot!(title="Composition estimations");
savefig(p, "composition-ages.svg");

# Write resampled to an HDF5 file
h5open("mcign.h5", "w") do file
    write(file, "mcign", mcign)
end

# Write model to an HDF5 file
h5open("pc1_model.h5", "w") do file
    write(file, "projection", projection(model))
end

## --- End of file
