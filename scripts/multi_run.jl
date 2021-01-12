"""
MPI script for running many inversions with random systematic biases applied to Crust1.0 data. 
Collects mean upper, middle, and lower compositions for each run and writes them to an output file.

Requries at least two MPI nodes, but if running on one machine for test use MPI options like: -np 2 --oversubscribe

Run like: 
mpiexec <mpi options> julia src/systematic_bias_test.jl <script options> 

"""

using MPI
MPI.Init()

using ArgParse
using DelimitedFiles
using StatGeochem
using HDF5
include("../src/config.jl")
include("../src/utilities.jl")
include("../src/inversionModel.jl")
include("../src/invertData.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Data folder"
        arg_type= String
        required = true
    "--num_runs", "-n"
    	help = "How many times to invert"
    	arg_type = Int
    	default = 4
    "--systematic"
    	help = "Apply systematic bias to Crust1.0 on each test?"
    	arg_type = Bool
    	default = false
    "--name"
    	help = "Name of run / output folder"
    	arg_type = String 
    	required= true
    "--data_source"
    	help = "Source for seismic data"
        arg_type = String
        range_tester = x -> (x in ["Shen","Crust1.0"])
        default = "Crust1.0"
    "--crack"
    	help = "Mean of upper crust total porosity."
    	arg_type = Float64  
    	default = 0.007
    	range_tester = x -> ((x >= 0) && (x < .5))
    "--fraction_crack"
    	help = "% of the pore space that is cracks (very low aspect ratio). 
    			Be careful with this, cracks have much larger effect on seismic props."
    	arg_type = Float64  
    	default = 0.05
    	range_tester = x -> ((x >= 0) && (x <= 1))
    "--cracked_samples"
    	help = "How many samples do we apply cracking to?"
    	arg_type = Float64  
    	default = 1.0
    	range_tester = x -> ((x >= 0) && (x <= 1))
end 
parsed_args = parse_args(ARGS, s)
N = parsed_args["num_runs"]

# MPI setup 
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
n_workers = MPI.Comm_size(comm)-1

"""
Run inversions, sending mean composition of each layer for each inversion to head. 
"""
function worker() 
	println("worker rank $(rank)")
	results = fill(-1.0,(3,10)) # layer, element 
	requested = fill(-1.0, (1)) # -1 if stop, else 1 
	while true 
		# Blocking send of data 
		MPI.Send(results, 0, rank, comm)

		# Blocking recv of next task
		MPI.Recv!(requested, 0, rank+10000, comm)

		# Check for kill signal (first index == -1)
		if requested[1] == -1
			break 
		end

		# Run test 
		# TODO 
		println("worker $(rank) test $(requested[1])")

		### SET UP CRACKING ### 
		outputPath = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])"
		crackFile = outputPath*"$(rank)_crack_profile.csv"
		ignFile = "data/"*parsed_args["data_prefix"]*"/bsr_ignmajors_1.csv"
		ign, header = readdlm(ignFile, ',', header=true)
		n = size(ign, 1)
		liquid_weights = [parsed_args["cracked_samples"]/2, parsed_args["cracked_samples"]/2, # dry, water 
			0, 1- parsed_args["cracked_samples"]] # magma, no cracking 
		crack_porosity = parsed_args["fraction_crack"] * parsed_args["crack"]
		pore_porosity = (1-parsed_args["fraction_crack"]) * parsed_args["crack"]
		
		if sum(liquid_weights) != 1 
			error("Sum of crack weights does not equal 1. Is your cracked_samples > 1?")
		end

		profiles = Array{Crack,1}([random_cracking(crack_porosity, pore_porosity, liquid_weights) for i in 1:n])
		write_profiles(profiles, crackFile)

		models = makeModels(parsed_args["data_prefix"], modelType=RangeModel) # see inversionModel.jl. returns a ModelCollection 
		
		age_model = EarthChemAge(10, 3)

		upperDat, (upperCrustbase, upperAge, upperLat, upperLong) = getAllSeismic(6, n=nOriginal(), 
			ageModel=age_model, latlong=true) # returns rho, vp, vpvs, tc1, age
		middleDat, (middleCrustbase, middleAge, middleLat, middleLong) = getAllSeismic(7, n=nOriginal(), 
			ageModel=age_model, latlong=true)
		lowerDat, (lowerCrustbase, lowerAge, lowerLat, lowerLong) = getAllSeismic(8, n=nOriginal(), 
			ageModel=age_model, latlong=true)

		# Result data per geotherm bin 
		results_upper, errors_upper = estimateComposition(models, UPPER, upperDat...)
		results_middle, errors_middle  = estimateComposition(models, MIDDLE, middleDat...)
		results_lower, errors_lower = estimateComposition(models, LOWER, lowerDat...)
		results_all = [results_upper, results_middle, results_lower]


		# set output of test to results 
		si_index = findfirst(isequal("SiO2"), RESAMPLED_ELEMENTS) # result 
		co_index = findfirst(isequal("CO2"), RESAMPLED_ELEMENTS)

		for (l, result) in enumerate(results_all)
			good = .~ isnan.(sum(result[:,si_index:co_index], dims=2))
			good = good[:]
	
			results[l, :] = mean(result[good,si_index:co_index], dims=1)
		end 

		# Remove crack file 
		rm(crackFile)
	end 
end 


"""
While still tasks left, recieve results from workers and send tasks to workers. 
(send tasks individually in case of differently efficient workers. )
"""
function head() 
	fileName = "data/$(parsed_args["data_prefix"])/$(parsed_args["name"])/multi_run.h5"
	if ispath(fileName)
		error("Output file at $fileName exists.")
	end
	mkpath("data/$(parsed_args["data_prefix"])/$(parsed_args["name"])")

	output = fill(-1.0,(3, 10, N)) # Collect data. 

	# Initial recvs for workers that haven't yet worked
	reqs = Vector{MPI.Request}(undef, n_workers)
	for i in 1:n_workers
		reqs[i] = MPI.Irecv!(Array{Float64,2}(undef,(3,10)), i,  i, comm)
	end

	for run_n in 1:N
		# Blocking wait for any worker 
		(worker_i, status) = MPI.Waitany!(reqs)
		
		# Blocking send of next task from ign to worker 
		MPI.Send(fill(1, (1)), worker_i, worker_i+10000, comm) # data, dst, id, comm 

		# Replace request for this worker with appropriate section of result array 
		reqs[worker_i] = MPI.Irecv!(view(output,:,:,run_n), worker_i, worker_i, comm)
	end

	# Wait for last worker results 
	MPI.Waitall!(reqs)

	# Send kill signal 
	kills = Vector{MPI.Request}(undef, n_workers)
	for i in 1:n_workers
		kills[i] = MPI.Isend(fill(-1.0,(1)), i, i+10000, comm)
	end
	MPI.Waitall!(kills)

	# Save output 
	h5write(fileName, "results", output)
end 

if rank == 0
	head() 
else 
	worker() 
end 







