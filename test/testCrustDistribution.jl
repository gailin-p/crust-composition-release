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

function testBinBoundaries() 
	bin1 = crustDistribution.binBoundaries(1)
	bin10 = crustDistribution.binBoundaries(10)
	if (length(bin1) != 2)
		throw(AssertionError("one-bin case broken"))
	end 
	if length(bin10) != 11
		throw(AssertionError("10 bin case broken"))
	end 
end 

testBinnedGeotherm()