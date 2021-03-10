using Test 

# Test that we can load data 
rm("data/crustDistribution.jld")

include("../src/crustDistribution.jl")

@testset "crustDistribution" begin 

# After removing file, should have been re-saved. 
@test isfile("data/crustDistribution.jld")

## crustDistribution should export depths, latitudes, longitudes 
@testset "Test initialized exports" begin 
	@test length(crustDistribution.all_lats) == length(crustDistribution.all_longs)
	@test length(crustDistribution.all_lats) == size(crustDistribution.depth, 1)
	@test size(crustDistribution.depth, 2) == 4 # geotherm, upper, middle, lower 
end 

@testset "Test latitude weighting" begin 
	@test crustDistribution.latLongWeight(-70) == crustDistribution.latLongWeight(70)
	@test crustDistribution.latLongWeight(-70) < crustDistribution.latLongWeight(-60)
	@test crustDistribution.latLongWeight(0) == 1
end

geotherm_lim = (50, 60)
@testset "Test fetching functions" for i in 1:100:1000
	# Right size return even as we chang inputs 
	@test size(crustDistribution.getCrustParams(i), 1) == i
	uncertain = crustDistribution.getCrustParams(i, uncertain=true)
	@test size(uncertain, 1) == i
	binned = crustDistribution.getCrustParams(geotherm_lim..., i, uncertain=false)
	@test size(binned, 1) == i

	## When uncertain, no values should be the same. 
	@test size(uncertain, 2) == 4 #geotherm, upper, middle, lower 
	@testset for j in 1:4 
		@test length(unique(uncertain[:,j])) == length(uncertain[:,j])
	end 

	## When geotherm bin, all isotherms fall within limits 
	@test all(binned[:,1] .>= geotherm_lim[1])
	@test all(binned[:,2] .<= geotherm_lim[2])
end

min_geotherm = minimum(crustDistribution.depth[:,1])
max_geotherm = maximum(crustDistribution.depth[:,1])
@testset "Test geotherm binning" for i in [1, 2, 15]
	bins = crustDistribution.binBoundaries(i)
	@test length(bins) == i+1
	@test bins[1] == min_geotherm
	@test bins[i+1] == max_geotherm
	if i > 1
		@test (bins[2] - bins[1]) == (bins[3] - bins[2])
	end
end 











end 