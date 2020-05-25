### Plot the relationship between seismic properties and composition 

using ArgParse
using HDF5
using MAT
using Statistics
using StatsBase
using Plots; gr();
include("../src/bin.jl")

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

layers = ["Upper", "Middle", "Lower"]

p = plot(size=(800,800));

for crust in layers

    rho = h5read(parsed_args["resampled"], "Calc_"*crust*"_Rho")
    sio2 = h5read(parsed_args["resampled"], "SiO2")

    test = map(!isnan, rho .+ sio2)

    xmin = percentile(rho[test], .5)
    xmax = percentile(rho[test],99.5)

    c, m, e, devs = bin(rho[test], sio2[test], xmin, xmax,
        length(sio2)/length(ign["SiO2"]), 20)

    plot!(p, c, m, yerror=e.*2, ylabel="Percent SiO2", xlabel="Rho (kg/m^3)",
         label=crust, markerstrokecolor=:auto);
end

savefig(p, "../output/SiO2-rho.pdf")
