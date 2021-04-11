include("run_mcmc.jl")

using ArgParse
using HDF5

s = ArgParseSettings()
@add_arg_table! s begin
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/scratch/gailin/" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        required = true
    "-n"
    	help = "Num loops"
    	arg_type = Int64
    	default = 10
end
parsed_args = parse_args(ARGS, s)

# remove output file 
out = "mcmc_out.h5"
rm(out, force=true)

# Data to invert (mean)
invert = [41.678471017560426 6.909417469309636 6.095079345979234 1.7252766096948702 2728.7077338784757; 41.678471017560426 20.193064518935472 6.413912813979953 1.7366797024951082 2802.8213184333244; 41.678471017560426 32.32255339791592 6.911113915442597 1.75810779522217 2933.306619695041]

runner = MCMCRunner(PerplexConfig(parsed_args["perplex"], parsed_args["scratch"]));

### RUN 
n = parsed_args["n"]
samples = fill(NaN, (3,n,10))
accepted = fill(false, (3,n))
ll = fill(NaN, (3,n))
for i in 1:3
    s, a, l = run_mcmc(n, runner, DataSample(invert[i,:]...))
    samples[i,:,:] .= s
    accepted[i,:] .= a
    ll[i,:] .= l
end 

h5open(out, "w") do file
    write(file, "samples", samples)  
    write(file, "accepted", accepted) 
    write(file, "loglike", ll) 
end

