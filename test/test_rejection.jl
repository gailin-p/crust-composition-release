using Test

include("../src/inversionModel.jl")
include("../src/rejectionModel.jl")
include("../src/config.jl")

@testset "Rejection model" begin

n = 2

ign = fill(-1.0,(n,length(PERPLEX_ELEMENTS)))
ign[:,1] = 1:n
ign[:,2] = [45, 25]
seismic = Array{Float64, 2}(undef,(n,4))
seismic[:,1] = 1:n

seismic[:,2] = [2000,1500] # 2 test populations: one within 10% of 1500, one within 10% of 2000
seismic[:,3] = [2, 1.5]
seismic[:,4] = [2, 1.5]  # 3 should be excluded from population 2 based on this value


model = RejectionModel(ign, seismic)


means, errs = estimateComposition(model, [2000., 1430., 500],
	[2., 1.6, .5],
	[2., 1.52, .5])

@test size(means) == (3,length(PERPLEX_ELEMENTS))
@test sum(isnan.(means[:,2])) == 0

println(means)

end
