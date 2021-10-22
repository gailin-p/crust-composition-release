using Test

include("../src/inversionModel.jl")
include("../src/rejectionModelSingle.jl")
include("../src/config.jl")

@testset "Single var rejection model" begin

n = 2

ign = fill(-1.0,(n,length(PERPLEX_ELEMENTS)))
ign[:,1] = 1:n
ign[:,2] = [45, 25]
seismic = Array{Float64, 2}(undef,(n,4))
seismic[:,1] = 1:n

seismic[:,2] = [2000,1500]
seismic[:,3] = [2, 1.5]
seismic[:,4] = [1, .8] # vpvs


vs_model = VsRejectionModel(ign, seismic)
@test(vs_model.seismic[:,2] == (seismic[:,3] ./ seismic[:,4]))

vp_model = VpRejectionModel(ign, seismic)
@test(vp_model.seismic[:,2] == seismic[:,3])

println(vs_model)



means, errs = estimateComposition(vs_model, [2000., 1430., 500],
	[2., 1.6, 2.1],
	[1., .75, .81])

@test size(means) == (3,length(PERPLEX_ELEMENTS))
@test sum(isnan.(means[:,2])) == 0

println(means)

end
