using MPI
MPI.Init()

using ArgParse
using DelimitedFiles

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for output data files"
        arg_type= String
        required = true
end 
parsed_args = parse_args(ARGS, s)

# MPI setup 
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
n_workers = MPI.Comm_size(comm)-1

# Global constants 
n = 3 # number of samples to send to each worker at once
sample_size = 16 # number of columns in ign per sample 


function worker() 
	println("worker rank $(rank)")
	results = fill(-1.0,(4,3,n)) # property (index, rho, vp, vpvs), layer, sample
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

		# Run perplex  TODO 
		results .= requested[1,1] # just for testing
	end
end 

"""
    head() 

Currently does not expect/allow worker failure 
"""
function head()
	println("head")
	# Load data used by head 
	ign = readdlm("data/"*parsed_args["data_prefix"]*"/bsr_ignmajors.csv", ',')
	ign = ign[1:20,:] # TODO just for testing
	n_samples = size(ign,1)
	output = fill(-1.0,(4,3,n_samples+n)) # Collect data. extra space for last worker run

	# Initial recvs for workers that haven't yet worked
	reqs = Vector{MPI.Request}(undef, n_workers)
	for i in 1:n_workers
		reqs[i] = MPI.Irecv!(Array{Float64,3}(undef,(4,3,n)), i,  i, comm)
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
		reqs[worker_i] = MPI.Irecv!(view(output,:,:,data_index:data_index+n), worker_i, worker_i, comm)
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
	print(output) # TODO temp 
end 


if rank == 0
	head() 
else 
	worker() 
end 
