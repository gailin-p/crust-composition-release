"""
Testing range type inversion model 
"""

using Statistics

include("../src/inversionModel.jl")

n = 5 

ign = Array{Float64, 2}(undef,(n,2))
ign[:,1] = 1:n
#ign[:,2] = 50 .+ 10*randn(n)
ign[:,2] = [45, 25, 30, 55, 45]
seismic = Array{Float64, 2}(undef,(n,4))
seismic[:,1] = 1:n
# Random test 
# seismic[:,2] = max.(1.0, 2500 .+ 500*randn(n))
# seismic[:,3] = max.(1.0, 4 .+ 2*randn(n))
# seismic[:,4] = max.(1.0, 2 .+ 2*randn(n))

seismic[:,2] = [2000,1500,1950,2100,2010] # 2 test populations: one within 10% of 1500, one within 10% of 2000
seismic[:,3] = [2, 1.5, 1.95, 2.1, 2.010]
seismic[:,4] = [2, 1.5, 0.5, 2.1, 2.010]  # 3 should be excluded from population 2 based on this value 


myModel = RangeModel(ign, seismic)
setError(myModel, (.1, .1, .1)) # 10% rel error 

# One with many samples, one with one sample, one too small, one too large 
means, errs = estimateComposition(myModel, [2000., 1000., 1600., 8000], [2.1, 1., 1.5, 8], [2., 1., 1.6, 8])

expected = [mean([45,55,45]), NaN, 25, NaN]
expected_std = [std([45,55,45]), NaN, NaN, NaN]

println("TEST RESULTS _____________")
println() 
println("means $means, errs $errs")
println("expected $expected, errs $expected_std")