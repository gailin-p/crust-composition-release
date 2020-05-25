using ArgParse
using StatsBase
using HDF5
using Plots; gr();

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Data folder"
        arg_type= String
        required = true
end 
parsed_args = parse_args(ARGS, s)

dat = h5read("data/$(parsed_args["data_prefix"])/results_systematic_bias_test.h5", "results")

upper = histogram(dat[1,1,:], xlims=(50,75), title="Upper", legend=false, grid=false, normalize=:pdf, titlefontsize=12, yaxis=false);
middle = histogram(dat[2,1,:], xlims=(50,75), title="Middle", legend=false, grid=false, normalize=:pdf, titlefontsize=12, yaxis=false);
lower = histogram(dat[3,1,:], xlims=(50,75), title="Lower", legend=false, grid=false, normalize=:pdf, titlefontsize=12,
	 yaxis=false, xlabel="% SiO2");

plot(upper, middle, lower, layout=(3,1), size=(400, 500));
savefig("data/$(parsed_args["data_prefix"])/output/systematic_bias_hist_layers.pdf")