################################
#
#  Averaging compositions before finding their seismic properties
#  is not the same as averaging those seismic properties afterwards.
#  How different, on average, are the two? What kind of error/uncertainty
#  might this introduce?
#
################################

using MPI
MPI.Init()

using ArgParse
using Plots; gr();
using JLD

include("../src/commutativityOfAverage.jl")

# Read user args: prefix for result files, input file names
s = ArgParseSettings()
@add_arg_table s begin
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/scratch/gailin/" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        required = true
    "--ign", "-i"
        help = "Path to ignmajors file"
        arg_type = String
        default = "ignmajors.csv"
    "--n_samples", "-n"
        help = "Number of samples to average"
        arg_type = Int
        default = 2
    "--n_runs", "-r"
        help = "How many sample pairs to run"
        arg_type = Int
        default = 2000
end
parsed_args = parse_args(ARGS, s)
r_requested = parsed_args["n_runs"]
n = parsed_args["n_samples"]

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nworkers = MPI.Comm_size(comm)-1


if rank != 0 # I'm not boss; run perplex then send my results to 0.
    # send
    # How many samples should I run?
    comm = MPI.COMM_WORLD
    r = ceil(Int, r_requested/nworkers)

    # data: average seismic properties, seismic properties of average, indices
    # for each of r runs, store upper, middle, lower; rho, vp, vpvs
    # dimensions: run, layer, property
    ave_properties = Array{Float64,3}(undef,(3,3,r))
    props_of_ave = Array{Float64,3}(undef,(3,3,r))
    indices = Array{Int,2}(undef,(n,r)) # store n indices for each of r runs, for reconstruction!

    # run samples
    commutativityOfAverage.runSamples!(props_of_ave, ave_properties, indices,
                r, n,
                parsed_args["ign"], parsed_args["perplex"], parsed_args["scratch"])

    # Send my data
    sreq1 = MPI.Isend(ave_properties, 0, rank, comm)
    sreq2 = MPI.Isend(props_of_ave, 0, rank+1000, comm)
    sreq3 = MPI.Isend(indices, 0, rank+2000, comm)
    stats = MPI.Waitall!([sreq1, sreq2, sreq3])
end

if rank == 0 # I'm the boss, wait for results
    rworker = ceil(Int, r_requested/nworkers) # runs per worker
    r = nworkers * rworker # total runs

    ave_properties = Array{Float64,3}(undef,(3,3,r))
    props_of_ave = Array{Float64,3}(undef,(3,3,r))
    indices = Array{Int,2}(undef,(n,r)) # store n indices for each of r runs, for reconstruction!

    # build recieve requests
    reqs = Vector{MPI.Request}()
    for w in 1:nworkers # worker ranks start at 1
        start = (w-1)*rworker+1
        rreq1 = MPI.Irecv!(view(ave_properties, :, :, start:start+rworker-1), w,  w, comm)
        rreq2 = MPI.Irecv!(view(props_of_ave, :, :, start:start+rworker-1), w,  w+1000, comm)
        rreq3 = MPI.Irecv!(view(indices, :, start:start+rworker-1), w,  w+2000, comm)
        push!(reqs, rreq1)
        push!(reqs, rreq2)
        push!(reqs, rreq3)
    end

    # Wait patiently
    MPI.Waitall!(reqs)

    # Save data and plots
    # Plot diffs
    diffs = props_of_ave-ave_properties
    layers = ["upper","middle","lower"]
    props = ["rho","vp","vpvs"]
    for l in 1:3
        for p in 1:3
           h = histogram(diffs[l,p,:], legend=false, xlabel="$(props[p]) difference")
           savefig(h, "output/average/diff_hist_$(props[p])_$(layers[l])_n$(n)_r$(r).pdf")
        end
    end

    # println(ave_properties)
    # println(props_of_ave)
    save("data/commutativity-n$(n)-r$(r).jld","ave_properties", ave_properties, "props_of_ave", props_of_ave, "indices", indices)

end
