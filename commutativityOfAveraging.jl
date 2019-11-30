################################
#
#  Averaging compositions before finding their seismic properties
#  is not the same as averaging those seismic properties afterwards.
#  How different, on average, are the two? What kind of error/uncertainty
#  might this introduce?
#
################################

using ArgParse
using DelimitedFiles
using StatGeochem
using ProgressMeter
using Statistics
using JLD
using Plots; gr();

# Average geotherm and layer depths
# General options
perplexdir = "/Users/gailin/dartmouth/crustal_structure/Perple_X/"
scratchdir = "/Users/gailin/dartmouth/crustal_structure/perplexed_pasta/"
exclude = ""
dataset = "hpha02ver.dat"
elts = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]
dpdz = 2900. * 9.8 / 1E5 * 1E3
depth = 40.2 # mean 550C depth of all samples
geotherm = 550.0/depth/dpdz
solutions_nof = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\n"
npoints = 20
upperStart = 1.08962179685088
upperBase = 13.442666718122878
middleBase = 25.551929607903677
lowerBase = 36.57912627354122
P_range = [1, ceil(Int,lowerBase*dpdz)] # only run perplex to bottom of crust for this geotherm

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

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nworkers = MPI.Comm_size(comm)-1

if rank != 0 # I'm not boss; run perplex then send my results to 0.
    # send
    # How many samples should I run?
    comm = MPI.COMM_WORLD
    r = ceil(Int, r_total/nworkers)

    # data: average seismic properties, seismic properties of average, indices
    # for each of r runs, store upper, middle, lower; rho, vp, vpvs
    # dimensions: run, layer, property
    ave_properties = Array{Float64,3}(undef,(r,3,3))
    props_of_ave = Array{Float64,3}(undef,(r,3,3))
    indices = Array{Int,2}(undef,(r,n)) # store n indices for each of r runs, for reconstruction!

    # run samples
    runSamples!(props_of_ave, ave_properties, indices,
                r, parsed_args["n_samples"],
                parsed_args["ign"], parsed_args["perplex"], parsed_args["scratch"])

    # Send my data
    sreq1 = MPI.Isend(ave_properties, 0, rank, comm)
    sreq2 = MPI.Isend(props_of_ave, 0, rank+1000, comm)
    sreq3 = MPI.Isend(indices, 0, rank+2000, comm)
    stats = MPI.Waitall!([sreq1, sreq2, sreq3])
end

if rank == 0 # I'm the boss, wait for results
    rworker = ceil(Int, r_total/nworkers) # runs per worker
    r = nworkers * rworker # total runs

    ave_properties = Array{Float64,3}(undef,(r,3,3))
    props_of_ave = Array{Float64,3}(undef,(r,3,3))
    indices = Array{Int,2}(undef,(r,n)) # store n indices for each of r runs, for reconstruction!

    # build recieve requests
    reqs = Vector{MPI.Request}()
    for w in 1:nworkers # worker ranks start at 1
        start = (i-1)*rworker+1
        rreq1 = MPI.Irecv!(view(ave_properties, start:start+rworker, :, :), w,  w, comm)
        rreq2 = MPI.Irecv!(view(props_of_ave, start:start+rworker, :, :), w,  w+1000, comm)
        rreq3 = MPI.Irecv!(view(indices, start:start+rworker, :), w,  w+2000, comm)
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
           h = histogram(diffs[:,l,p], legend=false, xlabel="$(props[p]) difference")
           savefig(h, "output/average/diff_hist_$(props[p])_$(layers[l])_n$(n)_r$(r).pdf")
        end
    end

    # println(ave_properties)
    # println(props_of_ave)
    save("data/commutativity-n$(n)-r$(r).jld","ave_properties", ave_properties, "props_of_ave", props_of_ave, "indices", indices)

end


### Run r samples
# Return (ave properties, properties of averages, indices)
# properties arrays are rx3x3, run x layer x property
function runSamples(r::Int, n::Int, ignFile::String, perplexdir::String, scratchdir::String)
    dat = readdlm(ignFile,',')

    # run r sample pairs through perplex
    @showprogress 1 "Running samples... " for i in 1:r
        allComp = Array{Float64,2}(undef,(n,10))
        allSeismic = Array{Any,1}(undef,n)

        # Run perplex for each composition
        for j in 1:n
            index = convert(Int, ceil(rand()*size(dat)[1]))
            indices[i,j] = index
            sample = dat[index,:]
            comp = sample[2:11]
            allComp[j,:] = comp
            # Run perple_x
            perplex_configure_geotherm(perplex, scratch, comp, elements=elts,
                P_range=P_range, geotherm=geotherm, dataset=dataset, solution_phases=solutions_nof,
                excludes="", npoints=npoints)

            seismic = perplex_query_seismic(perplexdir, scratchdir)
            allSeismic[j]=seismic
        end

        # Average seismic properties over samples (assume p is same at each index)
        # keys: "vp,km/s", "rho,kg/m3", "T(K)", "elements", "P(bar)", "vp/vs"
        aveVp = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
        aveVpVs = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))
        aveRho = Array{Float64, 1}(undef, length(allSeismic[1]["P(bar)"]))

        for p_index in 1:length(allSeismic[1]["P(bar)"])
            aveVp[p_index] = nanmean(map(y -> (allSeismic[y]["vp,km/s"][p_index]), 1:length(allSeismic)))
            aveVpVs[p_index] = nanmean(map(y -> (allSeismic[y]["vp/vs"][p_index]), 1:length(allSeismic)))
            aveRho[p_index] = nanmean(map(y -> (allSeismic[y]["rho,kg/m3"][p_index]), 1:length(allSeismic)))
        end

        # Run perplex on average composition
        aveComp = Array{Float64,1}(undef,10)
        for j in 1:10
            aveComp[j] = mean(allComp[:,j])
        end
        perplex_configure_geotherm(perplexdir, scratchdir, aveComp, elements=elts,
            P_range=P_range, geotherm=geotherm, dataset=dataset, solution_phases=solutions_nof,
            excludes="", npoints=npoints)

        seismicOfAve = perplex_query_seismic(perplexdir, scratchdir)

        # Average seismic properties over depth layer
        # (TODO is this ok? we're trying to *test* the legitimacy of averaging)
        upperTest = seismicOfAve["P(bar)"] .< upperBase * dpdz
        middleTest = (seismicOfAve["P(bar)"] .> upperBase * dpdz) .& (seismicOfAve["P(bar)"] .< middleBase * dpdz)
        lowerTest = (seismicOfAve["P(bar)"] .> middleBase * dpdz) .& (seismicOfAve["P(bar)"] .< lowerBase * dpdz)
        tests = [upperTest, middleTest, lowerTest]

        for l in 1:3 # upper, middle, lower
            ave_properties[i,l,1] = nanmean(aveRho[tests[l]])
            ave_properties[i,l,2] = nanmean(aveVp[tests[l]])
            ave_properties[i,l,3] = nanmean(aveVpVs[tests[l]])
            props_of_ave[i,l,1] = nanmean(seismicOfAve["rho,kg/m3"][tests[l]])
            props_of_ave[i,l,2] = nanmean(seismicOfAve["vp,km/s"][tests[l]])
            props_of_ave[i,l,3] = nanmean(seismicOfAve["vp/vs"][tests[l]])
        end
    end

    return (props_of_ave, ave_properties, indices)
end
