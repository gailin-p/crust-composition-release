include("../src/utilities.jl")
using Test 

@testset "Utilities tests" begin 

@testset "Test running mean" begin  
	runner = RunningMean()
	nums = randn(100)
	for n in nums 
		mean!(runner, n)
	end 
	@test abs(runner.m - mean(nums)) < 1e-16
end 

@testset "Test normalize" begin 
	a = [1 2 3; .1 .2 .3; .5 .5 .5]
	normalizeComp!(a)

	for row in 1:size(a,1)
		@test abs(sum(a[row,:]) - 100) < 1e-13
	end 
end


end 