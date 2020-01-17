include("../src/utilities.jl")

function testRunningMean() 
	runner = RunningMean()
	nums = randn(100)
	for n in nums 
		mean!(runner, n)
	end 
	if abs(runner.m - mean(nums)) > 1e-16
		throw(AssertionError("Running mean $(runner.m) does not equal mean $(mean(nums))"))
	else 
		println("Running mean test passed")
	end 
end 

testRunningMean()