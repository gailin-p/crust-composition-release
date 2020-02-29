using Statistics
using DelimitedFiles 


"""
Keep a running mean
m mean so far 
n number of samples so far 
"""
mutable struct RunningMean
	m::Float64
	n::Int
end 

function RunningMean()
	return RunningMean(0., 0)
end 

"""
Add a new number to a running mean 
"""
function mean!(running::RunningMean, new::Number)
	running.m = (running.n*running.m + new)/(running.n + 1)
	running.n += 1 
end 

"""
Normalize each row of input matrix (in place)
"""
function normalizeComp!(a::AbstractArray{Float64, 2})
	sums = sum(a,dims=2)
	for row in 1:size(a,1)
		a[row,:] .= (a[row,:] ./ sums[row]) .* 100
	end 
end 

"""
Write options to option file
"""
function writeOptions(filename, options)
	out = fill("", (length(keys(options)), 2))
	for (k, key) in enumerate(keys(options))
		out[k, 1] = key 
		out[k, 2] = string(options[key])
	end 
	writedlm(filename, out, ",")
end 
















