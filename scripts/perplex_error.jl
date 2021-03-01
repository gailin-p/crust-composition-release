"""
Compare perplex-calculated values for upper, middle, and lower crust 
with values from Hacker and Abers 2016. 
"""

using ArgParse
using DelimitedFiles
using StatGeochem
using JLD
include("../src/config.jl")
include("../src/seismic.jl")
include("../src/parsePerplex.jl")
include("../src/cracks.jl")
include("../src/crustDistribution.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/users/gailin/dartmouth/crustal_structure/perplexed_pasta/" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        default = "/users/gailin/resources/perplex-stable/"
    "--perplex_dataset"
    	help = "Which perplex thermodynamic dataset to use? eg hpha02ver.dat or hp02ver.dat"
    	arg_type = String
    	default = "hpha11ver.dat"
    "--solution_list"
    	help = "Use old or new solutions?"
    	arg_type = String
    	range_tester = x -> (x in ["old","new"]) 
    	default = "old"
end 
parsed_args = parse_args(ARGS, s)
perplex = parsed_args["perplex"]
scratch = parsed_args["scratch"]

fluid_endmembers = "abL\nanL\ndiL\nenL\nfaL\nfliq\nfoL\nkspL\nmliq\nqL\nsiL\nq8L\nfa8L\nfo8L\nsil8L\nh2oL\nh2o8L\n"

if parsed_args["solution_list"] == "old"
	# solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)"*
	# 		"\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar"*
	# 		"\nDo(HP)\n"
	solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)"*
	 		"\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar"*
	 		"\nDo(HP)\n"
	excludes = fluid_endmembers
else # From brenhin 
	solutions = "Augite(G)\nOpx(JH)\ncAmph(G)\noAmph(DP)\nfeldspar_B\nO(JH)\nSp(JH)\nGrt(JH)\nMica(W)\nBio(TCC)"*
			"\nChl(W)\nCtd(W)\nCrd(W)\nSa(WP)\nSt(W)\nIlm(WPH)\nAtg(PN)\nT\nB\nF\nDo(HP)\nScap\nChum\n"
	excludes = fluid_endmembers*"ged\nfanth\ngl\nilm\nilm_nol\n"
end

# "SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2"
# From kern dabie  
comps, h = readdlm("data/kern_dabie_comp.csv", ',', header=true)
h = h[:]
sample_names = comps[:,1][:]
comp_compat = zeros((size(comps,1),length(COMPOSITION_ELEMENTS)))
for (j, name) in enumerate(COMPOSITION_ELEMENTS)
    if name == "H2O_Total"
        name = "H2OC"
    end
    if name == "FeO"
        feoi = findfirst(isequal("FeO"),h)
        fe2o3i = findfirst(isequal("Fe2O3"), h)
        feo = [feoconversion(comps[i, feoi], comps[i, fe2o3i]) for i in 1:size(comps,1)]
        comp_compat[:,j] .= feo
    else
        comp_compat[:, j] .= comps[:,findfirst(isequal(name), h)]
    end
end



### Choose some conditions 
# median depth to 550 isotherm 
isotherm, depth = crustDistribution.getFormationParams()
println("Using formation params: depth to 550 $isotherm, depth $depth")
#isotherm = 39.0;
dpdz = 2900. * 9.8 / 1E5 * 1E3; # bar/km 
dtdz = 550.0/isotherm # K/km 
# perplex at 1/2 median depth of crust
#depth = 18
# calc t and p 
tts = depth * dtdz # linear from 0 at surface to 550 at isotherm.
pps = depth * dpdz 
geotherm = dtdz/dpdz
P_range = [depth*(9/10)*dpdz, depth*(11/10)*dpdz] # we only need a small range now 

# For seismic 
st = SeismicTransform()
# Temperature and pressure conditions and associated headers in Dabie data 
condition_map = Dict([("25 Mpa", (25, 20)), ("50 Mpa", (50, 20)), ("100 Mpa", (100, 20)), 
	("200 Mpa", (200, 20)), ("400 Mpa", (400, 20)), ("600 Mpa", (600, 20)),
	("100 C", (600, 100)), ("200 C", (600, 200)),  ("400 C", (600, 400)), ("600 C", (600, 600))])

out_header = ["sample" "P" "T" "porosity" "perplex rho" "perplex vp" "perplex vp/vs" "perplex mu" "perplex ks" "dabie rho" "dabie vp" "dabie vp/vs" "dabie mu" "dabie ks"]
out = fill(NaN, (size(comp_compat,1)*length(keys(condition_map)), length(out_header)))
out_endmembers = Dict() 

# Use dabie data as comparison 
dabie_vp, vp_header = readdlm("data/dabie_vp.csv", ',', header=true)
dabie_vs, vs_header = readdlm("data/dabie_vs.csv", ',', header=true)
dabie_rho, rho_header = readdlm("data/kern_dabie_seismic.csv", ',', header=true)
dabie_rho = dabie_rho[:,1:2] # first row is names, second row is density 


for i in 1:size(comp_compat,1)
	# Run perplex
	comp = comp_compat[i,:]
	perplex_configure_geotherm(perplex, scratch, comp, PERPLEX_COMPOSITION_ELTS,
	    P_range, 273.15, geotherm, dataset=parsed_args["perplex_dataset"], solution_phases=solutions,
	    excludes=fluid_endmembers, npoints=20)
	point = perplex_query_point(perplex, scratch, pps)
	properties = get_system_props(point)

	# Change temperature and pressure conditions of seismic properties 
	try 
		endmembers = parse_perplex_point(point)
		out_endmembers[i] = endmembers 
		j = 1
		for (k, v) in condition_map

			idx = (i-1)*length(keys(condition_map)) + j 
			if !all(isnan.(out[idx,:]))
				println("Indexing messed up")
			elseif (sample_names[i] != dabie_rho[i,1]) 
				println("Sample loc messed up")
			end 

			out[idx, 1] = i
			out[idx, 2:3] .= v
			
			# Seismic adjusted to this T, P 
			rho, vp, vs = [get_seismic(v[2]+273.15, v[1]*10, properties, endmembers, st)...]
			cracks = CrackProfile(v[1]*10/dpdz) # default cracking for this depth 
			out[idx, 4] = cracks.porosity
			out[idx,5:7] .=  apply_cracking(rho, vp, vp/vs, cracks) # returns rho, vp, vp/vs 
			out[idx,8:9] .= ks_mu(rho, vp, vs)
			
			# Find appropriate data from 3 different Dabie tables 
			out[idx, 10] = dabie_rho[i,2]
			d_vp = dabie_vp[(dabie_vp[:,1] .== sample_names[i]) .& (dabie_vp[:,2] .== "Mean"), 
					findfirst(isequal(k),vp_header[:])][1]
			d_vs = dabie_vs[(dabie_vs[:,1] .== sample_names[i]) .& (dabie_vs[:,2] .== "Mean"), 
					findfirst(isequal(k),vs_header[:])][1] 
			out[idx, 11] = d_vp 
			out[idx, 12] = d_vp/d_vs 
			out[idx, 13:14] .= ks_mu(dabie_rho[i,2]*1000, d_vp, d_vs)
			
			j += 1 
		end
	catch e
		if (isa(e, ParsePerplexError) | isa(e, SeismicError))
			println("\r\n\r\nCannot process point due to \r\n $e")
		else 
			throw(e)
		end 
	end
	#throw(RuntimeError()) # if want to die after 1 calc 
	#out[i,4:6] .= [properties["Density(kg/m3)"], properties["Vp(km/s)"], properties["Vs(km/s)"]]
end 

save("data/perplexEndmembers.jld", "dat", out_endmembers)

out = vcat(out_header, out)
writedlm("data/adjustedPerplexRes.csv", out, ',')

