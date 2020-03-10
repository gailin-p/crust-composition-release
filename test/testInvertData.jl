using Test

include("../src/invertData.jl")

function testEarthChemAgeModel()
	nbins = 10
	model = EarthChemAge(nbins, 1)

	# Test bins 
	ageBins = ageBins(model)
	@test length(ageBins) == nbins+1 # ageBins is top/bottom of bins 

	# Test age constraints 
end 