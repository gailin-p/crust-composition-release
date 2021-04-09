using Test 

include("../src/linearModel.jl")

@testset "Linear model" begin 

# index, rho, vp, vpvs 
seismic = [1. 10. 1^2 1.1
       	   2. 20. 2^2 1.
           3. 30. 3^2 1.1
           4. 40. 4^2 0.9
           5. 50. 5^2 1.2]

# 2 composition dependent vars 
ign = [(1. +1^2) 1. 
	(2. +2^2) 2. 
	(3. +3^2) 3. 
	(4. +4.5^2) 4.
	(5.5 +4^2) 5. ]

lm = LinearModel(ign, seismic)

@test size(lm.fit) == (2,4) # 2 dependent vars, 3 independent vars plus intercept 

@test abs(lm.fit[2,2]) < 1e-10 # second dependent var not dependent on x^2 component 

res, er = estimateComposition(lm, [25.,40.],[2.5,4.],[1.,1.])

@test size(res) == (2,2)

end 