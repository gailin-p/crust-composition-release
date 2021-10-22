"""
Compare perplex-calculated values for upper, middle, and lower crust
with values from Hacker and Abers 2016.
"""

using ArgParse
using DelimitedFiles
using StatGeochem
using JLD
using HDF5
include("../src/config.jl")
include("../src/seismic.jl")
include("../src/parsePerplex.jl")
include("../src/cracks.jl")
include("../src/crustDistribution.jl")
include("../src/utilities.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/Users/f0043n9/dartmouth/crustal_structure/perplex_scratch" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        default = "/Users/f0043n9/dartmouth/crustal_structure/perplex_687"
    "--perplex_dataset"
    	help = "Which perplex thermodynamic dataset to use? eg hpha02ver.dat or hp02ver.dat"
    	arg_type = String
    	default = "hpha11ver.dat"
    "--solution_list"
    	help = "Use old or new solutions?"
    	arg_type = String
    	range_tester = x -> (x in ["old","new"])
    	default = "old"
    "-n"
    	help = "How many T/P combos to try? if 1, use default 20 km / 450 C, else random"
    	arg_type = Int64
    	range_tester = x -> (x > 0)
    	default = 1
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


comp_compat, sample_names = convert_dabie("resources/dabie/kern_dabie_comp.csv")

# For seismic
st = SeismicTransform()
# Temperature and pressure conditions and associated headers in Dabie data
condition_map = Dict([("25 Mpa", (25, 20)), ("50 Mpa", (50, 20)), ("100 Mpa", (100, 20)),
	("200 Mpa", (200, 20)), ("400 Mpa", (400, 20)), ("600 Mpa", (600, 20)),
	("100 C", (600, 100)), ("200 C", (600, 200)),  ("400 C", (600, 400)), ("600 C", (600, 600))])


## Output if only 1 temp/pressure (csv)
nsamples = size(comp_compat,1)*length(keys(condition_map))
out_header = ["sample" "P" "T" "porosity" "perplex rho" "perplex vp" "perplex vp/vs" "perplex mu" "perplex ks" "dabie rho" "dabie vp" "dabie vp/vs" "dabie mu" "dabie ks"]
out_endmembers = Dict()
out = fill(NaN, (nsamples, length(out_header)))

## Output if many endmembers (hdf5)
many_out = fill(NaN, (parsed_args["n"], nsamples, 13))

# Use dabie data as comparison
dabie_vp, vp_header = readdlm("resources/dabie/dabie_vp.csv", ',', header=true)
dabie_vs, vs_header = readdlm("resources/dabie/dabie_vs.csv", ',', header=true)
dabie_rho, rho_header = readdlm("resources/dabie/kern_dabie_seismic.csv", ',', header=true)
dabie_rho = dabie_rho[:,1:2] # first row is names, second row is density

for n in 1:parsed_args["n"]

	# Choose conditions
	if parsed_args["n"] == 1
		depth, dtdz = crustDistribution.getFormationParams()
	else
		depth = rand()*50 # uniform 0 to 50 km
		temp = 200 + rand()*400 # uniform 200-600 C
		dtdz = temp / depth
	end
	tts = depth * dtdz # linear from 0 at surface to 550 at isotherm.
	pps = depth * dpdz
	geotherm = dtdz/dpdz
	P_range = [depth*(9/10)*dpdz, depth*(11/10)*dpdz] # we only need a small range now
	println("$n: Using formation params: temp $(dtdz*depth) deg C, depth $depth km")

	# Now run
        out .= NaN
	for i in 1:size(comp_compat,1)
		# Run perplex
		comp = comp_compat[i,:]
		perplex_configure_geotherm(perplex, scratch, comp, PERPLEX_COMPOSITION_ELTS,
		    P_range, 273.15, geotherm, dataset=parsed_args["perplex_dataset"], solution_phases=solutions,
		    excludes=fluid_endmembers, npoints=20)
		point = perplex_query_point(perplex, scratch, pps)

		# Change temperature and pressure conditions of seismic properties
		try
			properties = get_system_props(point)
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
			println(e)
			# if (isa(e, ParsePerplexError) | isa(e, SeismicError))
			# 	println("\r\n\r\nCannot process point due to \r\n $e")
			# else
			# 	throw(e)
			# end
		end
		# Save if many conditions
          		    many_out[n, :, 1] .= depth
		    many_out[n, :, 2] .= depth * dtdz
		    many_out[n, :, 3:9] .= out[:,1:7]
		    many_out[n, :, 10] .= (out[:,6] ./ out[:,7]) # vs
		    many_out[n, :, 11] .= (out[:,6] .- out[:,11]) # vp er
		    many_out[n, :, 12] .= (out[:,7] .- out[:,12]) # vpvs er
		    many_out[n, :, 13] .= (out[:,6] ./ out[:,7] .- (out[:,11] ./ out[:,12])) # vs er
             	end
end

if parsed_args["n"] == 1
	save("data/perplexEndmembers.jld", "dat", out_endmembers)

	out = vcat(out_header, out)
	writedlm("data/adjustedPerplexRes_687.csv", out, ',')
else
	h5write("data/perplexStability.h5", "res", many_out)
end
