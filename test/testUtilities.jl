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

function testNormalize()
	a = [1 2 3; .1 .2 .3; .5 .5 .5]
	normalizeComp!(a)

	for row in 1:size(a,1)
		if abs(sum(a[row,:]) - 100) > 1e-13
			println(sum(a[row,:]))
			throw(AssertionError("Row does not sum to 100"))
		end 
	end 
	println("Normalize test passed")
end

testRunningMean()

testNormalize()

