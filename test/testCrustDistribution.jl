include("../src/crustDistribution.jl")

# Test getting random samples in bin 
function testBinnedGeotherm()
	samples = crustDistribution.getCrustParams(40,44,100)
	geotherms = samples[:,1]
	if minimum(geotherms) <= 40
		throw(AssertionError("Minimum value less than requested"))
	end
	if maximum(geotherms) > 44
		throw(AssertionError("Maximum value less than requested"))
	end
end

testBinnedGeotherm()