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
using Statistics
using ProgressMeter: @showprogress

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-o"
        help = "Folder for output data files"
        arg_type= String
        required = true
    "--data", "-d"
        help = "Path to EarthChem mat file"
        arg_type = String
        default = "igncn1.mat"
    "--num_samples", "-n"
        help = "How many samples to produce"
        arg_type = Int
        default = 100000
    "--weight", "-w"
        help = "Weight according to silica content? Allowed: silica"
        arg_type = String
        range_tester = x -> (x in ["","silica"])
        default = ""
end
parsed_args = parse_args(ARGS, s)

# Read in mat file
ign = matread(parsed_args["data"])
ign["elements"] = string.(ign["elements"]) # Since .mat apparently can"t tell the difference between "B" and "b"

# List of elements we want to export (primarily element oxides)
elements = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2","tc1Crust"]

# Set uncertainties
for e in elements
    ign["err"][e] = ign[e] .* (ign["err2srel"][e] / 2) # default
    ign["err"][e][isnan.(ign[e])] .= nanstd(ign[e]) # missing data gets std TODO is this ok?
end

# Clean up data before resampling
# Start with all Fe as FeO TODO before or after resample? TODO what error to use for binned FeO?
ign["FeO"] = feoconversion.(ign["FeO"], ign["Fe2O3"], ign["FeOT"], ign["Fe2O3T"])

# Set undefined elements to average. TODO is this right?
for e in elements
    ign[e][isnan.(ign[e])] .= nanmean(ign[e])
end


# Resample
# TODO should I resample weighting for uniform SiO2?
if (parsed_args["weight"] == "silica") # weight according to silica content
    k = invweight(ign["SiO2"], (nanmaximum(ign["SiO2"])-nanminimum(ign["SiO2"]))/100)

    # Probability of keeping a given data point when sampling:
    # We want to select roughly one-fith of the full dataset in each re-sample,
    # which means an average resampling probability <p> of about 0.2
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    resampled = bsresample(ign, parsed_args["num_samples"], elements, p)
else
    resampled = bsresample(ign, parsed_args["num_samples"], elements)
end

# Write ignmajors

# Create 2d array of out element data to export
outtable = Array{Float64,2}(undef, length(resampled["SiO2"]), length(elements)+1)
for i = 1:length(elements)
    outtable[:,i+1] = resampled[elements[i]]
end
outtable[:,1] = Array(1:length(resampled["SiO2"]))

# Reject samples with suspicious anhydrous normalizations
anhydrousnorm = sum(outtable[:,2:9], dims=2)
t = (anhydrousnorm .< 101) .& (anhydrousnorm .> 90)

println("Discarded $(size(outtable,1) - sum(t)) because of suspicious anhydrous normalizations")

# Write accepted samples to file
dir = parsed_args["data_prefix"]
mkpath("data/"*dir) # make if does not exist

path = "data/"*dir*"/bsr_ignmajors.csv"
writedlm(path, round.(outtable[t[:],:], digits=5), ",")
