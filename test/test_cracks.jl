include("../src/cracks.jl")
using Test 

"""
SAMPLE CRACKS 
"""
## some porosity 
some_needles = CrackProfile(needle_properties, 0, 0, .01, 0, "d", "needle")
# more porosity 
more_needles = CrackProfile(needle_properties, 0, 0, .10, 0, "d", "needle")
# needles and water
some_needles_wet = CrackProfile(needle_properties, water_rho, water_K, .01, 0, "w", "needle")
## some porosity, some crack 
some_needles_some_crack = CrackProfile(needle_properties, 0, 0, .01, 0.001, "d", "needle")
## some crack 
some_crack = CrackProfile(needle_properties, 0, 0, 0, 0.001, "d", "needle")
## 0 porosity 
no_crack = CrackProfile(needle_properties, 0, 0, 0, 0, "d", "needle")
## Programmed no crack 
no_crack_type = NoCrack()

"""
SAMPLE ROCKS 
"""
r1 = [2600, 6, 1.75] # something felsic 
r2 = [3400, 7.5, 1.75] # something mafic 
r3 = [1500, 8, 3] # something wild 

"""
cracks.jl works in terms of rho, vp, vpvs, 
but for testing it's easier to predict the effect on vs 
vs = vp/(vp/vs)
Take an array of [rho, vp, vpvps] and return vs
"""
function vs(vel_array)
	return vel_array[2]/vel_array[3]
end

"""
TESTS 
"""

@testset "Cracking tests" begin 

## No cracking doesn't change rocks 
@testset "No cracking" begin
	@test all([apply_cracking(r1..., no_crack_type)...] .== r1)
	@test all([apply_cracking(r1..., no_crack)...] .== r1)
end

## After cracking, vp, vs, rho decrease 
@testset "Crack decreases vp and rho" for r in [r1, r2, r3]
	@testset "For all crack styles" for c in [some_needles, some_needles_wet, some_crack, some_needles_some_crack]
		cracked = [apply_cracking(r..., c)...]
		@test cracked[1] < r[1]
		@test cracked[2] < r[2]
		@test vs(cracked) < vs(r)
	end
end

## Additional cracking further decreases vp vs rho 
@testset "More cracking results in more decrease for vp and rho" for r in [r1, r2, r3]
	needles = [apply_cracking(r..., some_needles)...]
	cracks = [apply_cracking(r..., some_crack)...]
	both = [apply_cracking(r..., some_needles_some_crack)...]
	@test all(both[1:2] .< needles[1:2])
	@test all(both[1:2] .< cracks[1:2])
	@test vs(both) < vs(needles)
	@test vs(both) < vs(cracks)
end

## Additional cracking further decreases vp vs rho 
@testset "More porosity results in more decrease for vp and rho" for r in [r1, r2, r3]
	needles = [apply_cracking(r..., some_needles)...]
	more = [apply_cracking(r..., more_needles)...]
	@test all(more[1:2] .< needles[1:2])
	@test vs(more) < vs(needles)
end

### Fluid sub increases rho, decreases vs because vs = sqrt(mu/rho) and mu does not change, 
# has complicated effect on vp not tested here.
@testset "Water increases vp, rho" for r in [r1, r2, r3]
	needles = [apply_cracking(r..., some_needles)...]
	water = [apply_cracking(r..., some_needles_wet)...]
	@test needles[1] .< water[1]
	@test vs(needles) > vs(water) 
end


@testset "Izv conversions" begin 
	gs = [5, 10, 15, 30]
	as = izvestiya_conversion(gs)
	@test all(as[1:2] .> 1)
	@test as[3] == 1
	@test as[4] < 1
end




end