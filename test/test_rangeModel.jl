"""
Testing range type inversion model 
"""

using Statistics
using Test 

include("../src/inversionModel.jl")
include("../src/config.jl")

@testset "Range model" begin 

n = 5 

ign = fill(-1.0,(n,length(PERPLEX_ELEMENTS)))
ign[:,1] = 1:n
ign[:,2] = [45, 25, 30, 55, 45]
seismic = Array{Float64, 2}(undef,(n,4))
seismic[:,1] = 1:n

seismic[:,2] = [2000,1500,1950,2100,2010] # 2 test populations: one within 10% of 1500, one within 10% of 2000
seismic[:,3] = [2, 1.5, 1.95, 2.1, 2.010]
seismic[:,4] = [2, 1.5, 0.5, 2.1, 2.010]  # 3 should be excluded from population 2 based on this value 


myModel = RangeModel(ign, seismic)
setError(myModel, (.1, .1, .1)) # 10% rel error 

# One with many samples, one with one sample, one too small, one too large 
means, errs = estimateComposition(myModel, [2000., 1000., 1600., 8000], 
	[2.1, 1., 1.5, 8], 
	[2., 1., 1.6, 8])

expected = [mean([45,55,45]), NaN, 25, NaN]

@test means[1,2] in [45, 55, 45] # these three samples match 1st seismic
@test means[3, 2] == 25 # only one sample matches 3rd seismic 

# Now change to using mean 
setMean(myModel, true)
# same samples 
means, errs = estimateComposition(myModel, [2000., 1000., 1600., 8000], 
	[2.1, 1., 1.5, 8], 
	[2., 1., 1.6, 8])

@test isapprox(means[1,2], mean([45, 55, 45])) # these three samples match 1st seismic
@test means[3, 2] == 25 # only one sample matches 3rd seismic 














end 