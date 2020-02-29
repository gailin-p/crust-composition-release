"""
MPI script for running perplex jobs on a cluster. 
Requries at least two MPI nodes, but if running on one machine for test use MPI options like: -np 2 --oversubscribe

Run like: 
mpiexec <mpi options> julia src/runPerplex.jl <script options> 
"""

using MPI
MPI.Init()

using ArgParse
using DelimitedFiles
using StatGeochem
using HDF5
include("config.jl")
include("utilities.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for output data files"
        arg_type= String
        required = true
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/scratch/gailin/" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        required = true
    "--geotherm_bin", "-b"
    	help = "Number of perplex bin if using bins this run. (Just controls which sample csv we look for)"
        arg_type = Int
        default = 1
    "--perplex_dataset"
    	help = "Which perplex thermodynamic dataset to use? eg hpha02ver.dat or hp02ver.dat"
    	arg_type = String
    	default = "hpha02ver.dat"
end 
parsed_args = parse_args(ARGS, s)
perplex = parsed_args["perplex"]
scratch = parsed_args["scratch"]

# MPI setup 
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
n_workers = MPI.Comm_size(comm)-1

# Global constants 
n = 20 # number of samples to send to each worker at once
sample_size = length(PERPLEX_ELEMENTS) # number of columns in ign per sample 

# General perplex options 
exclude = ""
dataset = parsed_args["perplex_dataset"]
dpdz = 2900. * 9.8 / 1E5 * 1E3
# For now, use fluid and throw away results with NaN or 0 seismic properties. 
solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\nF\n"
npoints = 20

# Perplex labels used by StatGeochem perplex interface 
prop_labels = ["rho,kg/m3","vp,km/s","vp/vs"]
p_label = "P(bar)"


function worker() 
	println("worker rank $(rank)")
	results = fill(-1.0,(2,4,3,n)) # result|error, property (index, rho, vp, vpvs), layer, sample
	requested = fill(-1.0, (n, sample_size))
	while true 
		# Blocking send of data 
		MPI.Send(results, 0, rank, comm)

		# Blocking recv of next data 
		MPI.Recv!(requested, 0, rank+10000, comm)
		println("Worker $(rank) got indices $(requested[:,1])")

		# Check for kill signal (first index == -1)
		if requested[1,1] == -1
			break 
		end

		# Run perplex  
		for i in 1:n 
			index = requested[i,1]
			comp = requested[i,2:length(COMPOSITION_ELEMENTS)+1]
			tc1 = requested[i,findfirst(isequal("geotherm"),PERPLEX_ELEMENTS)]
			layers = requested[i,findfirst(isequal("upper"),PERPLEX_ELEMENTS):findfirst(isequal("lower"),PERPLEX_ELEMENTS)]
			if index == -1 # out of real samples 
				break
			end

			# Sample-specific perplex options 
			geotherm = 550.0/tc1/dpdz
        	P_range = [1, ceil(Int,layers[3]*dpdz)] # run to base of crust. 

        	# Run perplex
        	perplex_configure_geotherm(perplex, scratch, comp, PERPLEX_COMPOSITION_ELTS,
                P_range, 273.15, geotherm, dataset=dataset, solution_phases=solutions,
                excludes="", index=rank, npoints=npoints)
            seismic = perplex_query_seismic(perplex, scratch, index=rank)

            # discard below first NaN or 0 value in any property 
            # find first NaN or 0 value across props
            bad = Inf 
			for prop_i in 1:length(prop_labels) 
            	prop = prop_labels[prop_i]
   				this_bad = findfirst(x->(isnan(x) | (x < 1e-6)), seismic[prop])
   				if !isnothing(this_bad)
   					bad = Int(min(this_bad, bad))
   				end
   			end
   			# Remove all values after first bad value 
   			if bad <= length(seismic[p_label]) # some values to be handled 
   				for prop in prop_labels # find first NaN or 0 value across props
					seismic[prop][bad:end] .= NaN 
				end
			end

            # Find per-layer mean for each property 
            p_layers = [l*dpdz for l in layers] # convert depths to pressures 
            pressure = seismic[p_label]
            upper = pressure .<= p_layers[1]
            middle = (pressure .> p_layers[1]) .& (pressure .<= p_layers[2])
            lower = (pressure .> p_layers[2]) .& (pressure .<= p_layers[3])
            for prop_i in 1:length(prop_labels)
            	prop = prop_labels[prop_i]
            	results[1,prop_i+1,1,i] = 1/(nanmean(1 ./ seismic[prop][upper]))
            	results[1,prop_i+1,2,i] = 1/(nanmean(1 ./ seismic[prop][middle]))
            	results[1,prop_i+1,3,i] = 1/(nanmean(1 ./ seismic[prop][lower]))
            	# Error 
            	results[2,prop_i+1,1,i] = nanstd(seismic[prop][upper])
            	results[2,prop_i+1,2,i] = nanstd(seismic[prop][middle])
            	results[2,prop_i+1,3,i] = nanstd(seismic[prop][lower])
            end

            # Set index in results 
            results[:,1,:,i] .= index

        end
	end
end 

"""
    head() 

Currently does not expect/allow worker failure 
"""
function head()
	println("head")
	writeOptions("data/"*parsed_args["data_prefix"]*"/runPerplex_options_$(parsed_args["geotherm_bin"]).csv", parsed_args)
	# Load data used by head 
	fileName = "data/"*parsed_args["data_prefix"]*"/bsr_ignmajors_$(parsed_args["geotherm_bin"]).csv"
	ign = readdlm(fileName, ',')
	n_samples = size(ign,1)
	output = fill(-1.0,(2,4,3,n_samples+n)) # Collect data. extra space for last worker run

	# Initial recvs for workers that haven't yet worked
	reqs = Vector{MPI.Request}(undef, n_workers)
	for i in 1:n_workers
		reqs[i] = MPI.Irecv!(Array{Float64,4}(undef,(2,4,3,n)), i,  i, comm)
	end

	for data_index in 1:n:n_samples
		# Blocking wait for any worker 
		(worker_i, status) = MPI.Waitany!(reqs)

		# Next n samples for this worker 
		to_send = fill(-1.0, (n, sample_size))
		if n_samples >= data_index + n 
			to_send = ign[data_index:data_index+n-1, :]
		else # running out, half-fill to_send
			remaining = n_samples - data_index + 1
			to_send[1:remaining,:] = ign[data_index:end,:]
		end 
		
		# Blocking send of next data from ign to worker 
		MPI.Send(to_send, worker_i, worker_i+10000, comm) # data, dst, id, comm 

		# Replace request for this worker with appropriate section of result array 
		reqs[worker_i] = MPI.Irecv!(view(output,:,:,:,data_index:data_index+n), worker_i, worker_i, comm)
	end

	# Wait for last worker results 
	MPI.Waitall!(reqs)

	# Send kill signal 
	kills = Vector{MPI.Request}(undef, n_workers)
	for i in 1:n_workers
		kills[i] = MPI.Isend(fill(-1.0,(n, sample_size)), i, i+10000, comm)
	end
	MPI.Waitall!(kills)

	# Save output 
	output = output[:,:,:,1:n_samples] # discard any trailing -1 
	fileName = "data/"*parsed_args["data_prefix"]*"/perplex_out_$(parsed_args["geotherm_bin"]).h5"
	h5write(fileName, "results", output[1,:,:,:])
	h5write(fileName, "errors", output[2,:,:,:])
end 


if rank == 0
	head() 
else 
	worker() 
end 
