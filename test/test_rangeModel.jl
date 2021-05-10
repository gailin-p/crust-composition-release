"""
Testing range type inversion model 
"""

using Statistics
using Test 

include("../src/inversionModel.jl")
include("../src/rangeModel.jl")
include("../src/config.jl")

@testset "Range model" begin 

n = 5

ign = fill(-1.0,(n,length(PERPLEX_ELEMENTS)))
ign[:,1] = 1:n
ign[:,2] = [45, 25, 30, 55, 46]
seismic = Array{Float64, 2}(undef,(n,4))
seismic[:,1] = 1:n

seismic[:,2] = [2000,1500,1950,2100,2010] # 2 test populations: one within 10% of 1500, one within 10% of 2000
seismic[:,3] = [2, 1.5, 1.95, 2.1, 2.010]
seismic[:,4] = [2, 1.5, 0.5, 2.1, 2.010]  # 3 should be excluded from population 2 based on this value 

# use small uncertainty so tests pass 
small = MvNormal(zeros(3), 1e-10 .* [1 0 0 ; 0 1 0; 0 0 1])

myModel = RangeModel(ign, seismic, small, relerr=.1) # 10% rel error 
#setError(myModel, (.1, .1, .1)) # 10% rel error 

# One with many samples, one with one sample, one too small, one too large, one that relies on % 
means, errs = estimateComposition(myModel, [2000., 1000., 1600., 8000, 2004.], 
	[2., 1., 1.5, 8, 2.01], 
	[2., 1., 1.6, 8, 2.01])

@test means[1,2] == 45 # the best of the three samples (45, 55, 46) match 1st seismic
@test means[3, 2] == 25 # only one sample matches 3rd seismic 
@test means[5, 2] == 46 # requires % distance eval to pass, scalar distance gives 45 as closest.
@test isequal(means[2, 2], NaN) # none matched 
@test isequal(means[4, 2], NaN) # none matched 

########### Using only vp 
small = MvNormal([0], [1e-10])
vpModel = RangeModel(ign, hcat(seismic[:,1], seismic[:,3]), small)
means, errs = estimateComposition(vpModel, [2., 1., 1.5, 8, 2.01])

@test means[1,2] == 45
@test isequal(means[2,2], NaN) # sample is too far 
@test means[3,2] == 25
@test isequal(means[4,2], NaN) 
@test means[5,2] == 46

######### Using Vp model builder 
vpModel2 = VpRangeModel(ign, seismic)
# Should look same 
@test vpModel.perms == vpModel2.perms
@test vpModel.lookups == vpModel2.lookups 


# # Now change to using mean 
# setMean(myModel, true)
# # same samples 
# means, errs = estimateComposition(myModel, [2000., 1000., 1600., 8000], 
# 	[2.1, 1., 1.5, 8], 
# 	[2., 1., 1.6, 8])

#@test isapprox(means[1,2], mean([45, 55, 46])) # these three samples match 1st seismic
#@test means[3, 2] == 25 # only one sample matches 3rd seismic 














end 