using Statistics


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