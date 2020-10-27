"""
Functions to calculate the seismic properties of rocks with cracks and fluid. 
Based on: 

David, E.C., Zimmerman, R.W., 2011. Elastic moduli of solids containing spheroidal pores. 
	Int J Eng Sci 49, 544–560. https://doi.org/10.1016/j.ijengsci.2011.02.001


Testing: 
jupyter notebook in visualization/porosity demonstrates that these replicate figures in David & Zimmerman

"""

using Random 
using StatsBase
using DelimitedFiles
using Distributions
using Roots


############### CRACKING BY DEPTH FUNCTION 

function kFromDepth(z) ### permeability from depth fn, from Manning & Ingebritsen 1999
	k = exp(-3.2*log(z)-14)
	return k 
end 

function porosityFromK(k)  #### porosity is some multiple of k 
	# for plutonic rocks, factors range from 10^2 to 10^15. Conservative: 10^2 to 10^4
	factor = rand()*(4-2) + 2
	return min(.015, k*10^factor)
end 

function porosityFromDepth(z)
	return porosityFromK(kFromDepth(z))
end 

############### END-TO-END FUNCTIONS 

abstract type Crack end 

# Each crack profile has a porosity for either needles or spheres, plus a crack porosity. 
struct CrackProfile <: Crack 
	cracking_fn::Function 
	fluid_density::Number # 0 for dry cracks 
	K_fluid::Number # 0 for dry cracks
	porosity::Number
	crack_porosity::Number
	liquid::String
	shape::String
end 

struct NoCrack <: Crack end 

"""
	props(CrackProfile)
Return base type fields of a crack profile for saving to a .csv
"""
function props(profile::CrackProfile)
	return [profile.fluid_density, profile.K_fluid, profile.porosity, profile.crack_porosity, profile.liquid, profile.shape]
end 

function props(profile::NoCrack)
	return [NaN, NaN, NaN, NaN, "none", "none"]
end 

# Load data 
# from Rivers and Carmichael 1987, "Ultrasonic Studies of Silicate Melts"
magma_dat, header = isfile("data/magma_moduli.csv") ? readdlm("data/magma_moduli.csv", ',', header=true) : readdlm("../data/magma_moduli.csv", ',', header=true)
magma_densities = magma_dat[:,3] .* 1000
magma_K = magma_dat[:, 5]
relerr_magma_K = .03 
relerr_magma_rho = .02

water_K = 2.8 # GPa
water_rho = 1000 # kg/m^3
relerr_water_K = .2 # Big uncertainties to allow for dissolved mierals 
relerr_water_rho = .2

"""
	random_cracking()
Return a CrackProfile for a random crack/fluid profile or no cracks. 
"""
function random_cracking(mean_crack_porosity::Number=.007, mean_pore_porosity::Number=.15,
	liquid_weights::Array{Float64, 1}=[.2, .2, .2, .4], # dry, water, magma, none 
	shape_weights::Array{Float64, 1}=[.5, .5]) # probability of sphere vs needle 
	# Dry, water or magma? 
	liquid_type = sample(["d","w","m","none"], Weights(liquid_weights))

	# select properties of liquid 
	if liquid_type == "d"
		K = 0
		rho = 0
	elseif liquid_type == "w"
		d_K = Normal(water_K, water_K*relerr_water_K)
		d_rho = Normal(water_rho, water_rho*relerr_water_rho)
		K = max(0, rand(d_K))
		rho = max(0, rand(d_rho))
	elseif liquid_type == "m"
		chosen_K = rand(magma_K)
		chosen_rho = rand(magma_densities)
		d_K = Normal(chosen_K, chosen_K*relerr_magma_K)
		d_rho = Normal(chosen_rho, chosen_rho*relerr_magma_rho)
		K = max(0, rand(d_K))
		rho = max(0, rand(d_rho))
	elseif liquid_type == "none"
		return NoCrack()
	end 

	# What shape pores? (name and function)
	cracking_fn = sample([(needle_properties, "needle"), (sphere_properties, "sphere")], Weights(shape_weights))


	# What porosity? 
	#porosity = rand() * .15 # uniform distribution of porosities, 0 to .15
	d_pore = LogNormal(log(mean_pore_porosity), 1)
	porosity = min(.5, rand(d_pore))

	# What cracking porosity? 
	d_porosity = LogNormal(log(mean_crack_porosity), 1) # Normal(mean_crack_porosity, .0002)
	crack_porosity = min(.1, rand(d_porosity))

	return CrackProfile(cracking_fn[1], rho, K, porosity, crack_porosity, liquid_type, cracking_fn[2])
end 

""""
Write and read profiles from a csv file -- 
allows for saving and later retrieving 
cracks for a model run. 
"""
function get_profiles(filename::String)
	cracking_functions = Dict([("needle", needle_properties), ("sphere", sphere_properties)])

	(dat, header) = readdlm(filename, ',', header=true)
	profiles = Array{Crack, 1}(undef, size(dat,1))
	for i in 1:size(dat,1)
		if dat[i,6] == "none"
			profiles[i] = NoCrack()
		else 
			profiles[i] = CrackProfile(cracking_functions[dat[i,6]], dat[i,:]...)
		end 
	end 
	return profiles 
end 

function write_profiles(profiles::Array{Crack, 1}, filename::String)
	header = string.(fieldnames(CrackProfile))[2:end] # Don't write first field, which is a function.
	output = Array{Union{String, Number}, 2}(undef, length(profiles)+1, length(header))
	output[1,:] = [header...]

	for i in 1:length(profiles)
		output[i+1,:] = props(profiles[i])
	end 
   
	writedlm(filename, output, ',')
end 

function apply_cracking(rho, vp, vpvs, profile::NoCrack)
	return rho, vp, vpvs
end 

""""
Apply cracking to these specific props 
"""
function apply_cracking(K0, mu0, rho0, fn::Function, porosity, fluid_density, fluid_K)
	poissonsd, KdK0 = fn(porosity, poissons(K0, mu0))
	Kd = KdK0*K0 # Dry k 
	mud = mu_from_poissons(Kd, poissonsd) # Dry mu 
	# Apply fluid (if K_fluid and fluid density are 0, this has no effect --> dry cracks)
	Kw = fluid_sub(K0, fluid_K, Kd, porosity)
	rho_w = wet_density(rho0, fluid_density, porosity)
	# Return input props but modified 
	return Kw, mud, rho_w
end 

function apply_cracking(rho, vp, vpvs, profile::CrackProfile)
	vs = vp/vpvs
	K0, mu0 = moduli_from_speeds(vp, vs, rho)

	if profile.porosity == 0 
		rho_w_pore = rho 
		vpw_pore = vp
		vsw_pore = vp / vpvs
	else 
	# try	
		### Effect on seismic velocity from pores (speheres and needles)
		Kw, mud, rho_w_pore = apply_cracking(K0, mu0, rho, 
			profile.cracking_fn, profile.porosity, 
			profile.fluid_density, profile.K_fluid)
		vpw_pore, vsw_pore = speed_from_moduli(Kw, mud, rho_w_pore)
	end 

	if profile.crack_porosity == 0 
		rho_w_crack = rho
		vpw_crack = vp 
		vsw_crack = vp / vpvs
	else 
		# # Repeat for cracking 
		Kw, mud, rho_w_crack = apply_cracking(K0, mu0, rho, 
			cracked_properties, profile.crack_porosity, 
			profile.fluid_density, profile.K_fluid)
		vpw_crack, vsw_crack = speed_from_moduli(Kw, mud, rho_w_crack)
	end 

	# velocity and density averages
	rho_w = 2/(1/rho_w_crack + 1/rho_w_pore)
	vpw = 2/(1/vpw_crack + 1/vpw_pore)
	vsw = 2/(1/vsw_crack + 1/vsw_pore)

	return rho_w, vpw, vpw/vsw
end 

"""
Modify a directory of perplex results by applying cracking to all. 
Use given list of cracking profiles. 
"""
function apply_cracking!(perplex_results::Array{Float64,3}, profiles::Array{Crack, 1}, upper::Bool) # output by perplex, dims property (index, rho, vp, vpvs), layer, sample 
	println("Applying cracking to all layers...")

	if length(profiles) != size(perplex_results, 3)
		error("Must be same number of cracking profiles as samples.")
	end 

	if upper
		println("applying cracking to upper crust only ")
		layers = [1]
	else 
		println("applying cracking to all layers")
		layers = 1:3
	end 

	for layer in layers
		println("layer $(layer)")
		for i in 1:size(perplex_results,3)
			if isnan(sum(perplex_results[2:4,layer,i]))
            	continue
        	end
			perplex_results[2:4,layer,i] = [apply_cracking(perplex_results[2:4,layer,i]..., profiles[i])...]
		end 
	end 
	return profiles
end 

############### CRACKS

CRACK_ALPHA = .01
MAX_CRACK_POROSITY = 4/3*pi*CRACK_ALPHA

"""
	gamma(porosity, alpha)
Crack density for a given porosity and crack aspect ratio.
"""
function gamma(porosity, alpha)
    return porosity / (4/3 * pi * alpha)
end

"""
	cracked_properties(porosity, solid poisson ratio)
Get new poisson's ratio and ratio of new/solid bulk moduli for cracks (assumed aspect ratio .01) 
"""
function cracked_properties(porosity, nu0)
    crack_density = gamma(porosity, CRACK_ALPHA)
    nu = nu0 * exp(-8/5 * crack_density)
    KK0 = ((1-2*nu0) * exp(-16/9 * crack_density)) / (1 - 2*nu0*exp(-8/5 * crack_density))
    #println("cracked $nu $KK0 from original $nu0")
    return nu, KK0
end

########## NEEDLES 

"""
	needle_properties(porosity, solid poisson ratio)
limit alpha = inf
""" 
function needle_properties(porosity, nu0)
    # Constants 
    lambda = (7-sqrt(29))/8
    n1 = 15/98 * (1+15/sqrt(29))
    N1 = 2/49 * (16 + 93/sqrt(29))
    
    nu = lambda - (lambda - nu0)*(1-porosity)^(1/n1)
    KK0 = (1-porosity)^(N1/n1) / ( (1-2*(lambda - (lambda-nu0)*(1-porosity)^(1/n1)))/(1-2*nu0))
    return nu, KK0
end 

########## SPHERES 

"""
Helper function for sphere_properties
"""
function build_nu_solver(porosity, nu0)
    return function nu_solver(nu)
        left = 1 - porosity 
        right1 = ((1-nu)/(1-nu0))
        right2 = ((1-5*nu)/(1-5*nu0))
        right3 = ((1+nu)/(1+nu0))
        if (right1 < 0) | (right2 < 0) | (right3 < 0)
            return 0.0
        end
        return right1^(-1/6) * right2^(5/6) * right3^(-2/3) - left
    end 
end 

"""
	sphere_properties(porosity, solid poisson ratio)
alpha = 1
""" 
function sphere_properties(porosity, nu0)
    nu = find_zero(build_nu_solver(porosity, nu0), nu0)
    term1 = Complex((1-5*nu)/(1-5*nu0))
    term2 = Complex((1-2*nu)/(1-2*nu0))
    term3 = Complex((1+nu)/(1+nu0))
    KK0 = real(term1^(5/3)) * real(term2^(-1)) * real(term3^(-2/3))
    return nu, KK0
end


############# FLUID 

"""
    fluid_sub(K0, Kf, Kd)
Gassmann's equations; formulation from Han and Batzle eq 2, 3
K0 = bulk moduli of mineral grain
Kf = bulk moduli of fluid
Kd = bulk moduli of dry frame 
returns Ks, bulk moduli of saturated frame (shear modulus doesn't change)
"""
function fluid_sub(K0, Kf, Kd, porosity)
    top = K0 * (1 - Kd/K0)^2
    bottom = 1 - porosity - Kd/K0 + porosity*(K0/Kf)
    return Kd + top/bottom
end 

function wet_density(rho0, rhof, porosity)
	return (1-porosity)*rho0 + porosity*rhof
end

"""
return K, mu (bulk modulus, shear modulus)
"""
function moduli_from_speeds(vp, vs, density)
    mu = density*vs^2
    K = density*vp^2 - (4/3)*mu
    return K, mu
end    

"""

"""
function speed_from_moduli(K, mu, density)
    vs = sqrt(mu/density)
    vp = sqrt((K + (4/3)*mu)/density)
    return vp, vs
end

function poissons(K, mu)
    return (3*K - 2*mu)/(6*K+2*mu)
end

function mu_from_poissons(K, poissons)
    top = poissons*6*K - 3*K
    bottom = -2 - 2*poissons
    return top/bottom
end



