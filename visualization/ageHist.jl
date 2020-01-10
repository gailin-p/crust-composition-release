### Compare age distribution of data before and after resampling.

using ArgParse
using MAT
using HDF5
using Plots; gr();

s = ArgParseSettings()
@add_arg_table s begin
    "--original", "-o"
        help = "Path to ign"
        arg_type = String
        required = true
    "--resampled", "-r"
        help = "Path to mcign"
        arg_type= String
        required = true
end
parsed_args = parse_args(ARGS, s)

ign = matread(parsed_args["original"])
ign_ages = ign["Age"]

mcign_ages = h5read(parsed_args["resampled"], "Age")

p1 = histogram(ign_ages, xlims=(0,4500), bins=60, xlabel="Age (Ma)", yaxis=nothing, title="Before resampling", legend=false, fillalpha=0)
p2 = histogram(mcign_ages, xlims=(0,4500), bins=60, xlabel="Age (Ma)", yaxis=nothing, title="After resampling", legend=false)

p = plot(p1, p2, layout=(1,2), size=(1000,700))
savefig(p, "../output/ageHist.pdf")
