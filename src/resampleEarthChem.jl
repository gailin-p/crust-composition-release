"""
Resample earthchem data set (original in igncn1.mat) to a perplex-ready ignmajors.csv

Advantages of this over using perplex directly on EarthChem samples:
- most EarthChem samples are missing H2O data
- some EarthChem samples behave badly at high pressure (0 vp) because Perplex
  predicts liquid, but discarding these samples
  at depth might bias the compositions we allow in the inversion for deeper layers.
  If all samples are resampled, we can discard at depth those that behave unphysically

Use:
Run from base folder (runPerplexBatchVp)
Places resampled file at data/<data_prefix>/ignmajors.csv
"""

using ArgParse
using StatGeochem
using MAT
using DelimitedFiles
using Logging
using ProgressMeter: @showprogress

s = ArgParseSettings()
@add_arg_table s begin
    "--data", "-d"
        help = "Path to EarthChem mat file"
        arg_type = String
        default = "igncn1.mat"
    "--data_prefix", "-o"
        help = "Folder for output data files"
        arg_type= String
        required = true
    "--num_samples", "-n"
        help = "How many samples to produce"
        arg_type = Int
        default = 100000
end
parsed_args = parse_args(ARGS, s)

# Read in mat file
ign = matread(parsed_args["data"])
ign["elements"] = string.(ign["elements"]) # Since .mat apparently can"t tell the difference between "B" and "b"

# List of elements we want to export (primarily out element oxides)
elements = ["index","SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2","tc1Crust"]

# Use default 2-sigma relative uncertainties from err2srel
for e in elements
    ign["err"][e] = ign[e] .* (ign["err2srel"][e] / 2)
end
# Use larger uncertainty for missing H2O and CO2
for e in ["H2O_Total", "CO2"]
    ign["err"][e][isnan.(ign[e])] .= nanstd(ign[e])
end

# Clean up data before resampling
# Start with all Fe as FeO TODO before or after resample? TODO what error to use for binned FeO?
ign["FeO"] = feoconversion.(ign["FeO"], ign["Fe2O3"], ign["FeOT"], ign["Fe2O3T"])

# Set undefined H2O and CO2 to average. TODO is this right?
ign["H2O_Total"][isnan.(ign["H2O_Total"])] .= nanmean(ign["H2O_Total"])
ign["CO2"][isnan.(ign["CO2"])] .= nanmean(ign["CO2"])


# Resample
# TODO should I resample weighting for uniform SiO2? 
resampled = bsresample(ign, parsed_args["num_samples"], elements)


# Write ignmajors

# Create 2d array of out element data to export
outtable = Array{Float64,2}(undef, length(resampled["SiO2"]), length(elements))
for i = 1:length(elements)
    outtable[:,i] = resampled[elements[i]]
end

# Reject samples with missing data # TODO shouldn't be any now I think...
t = .~ any(isnan.(outtable), dims=2)

# Reject samples with suspicious anhydrous normalizations
anhydrousnorm = sum(outtable[:,2:9], dims=2)
t = t .& (anhydrousnorm .< 101) .& (anhydrousnorm .> 90)

println("Discarded $(size(outtable,1) - size(t,1))")

# Write accepted samples to file
dir = parsed_args["data_prefix"]
mkpath("data/"*dir) # make if does not exist

path = "data/"*dir*"bsr_ignmajors.csv"
writedlm(path, round.(outtable[t[:],:], digits=5), ",")
