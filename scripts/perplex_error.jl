"""
Compare perplex-calculated values for upper, middle, and lower crust 
with values from Hacker and Abers 2016. 
"""

using ArgParse
using DelimitedFiles
using StatGeochem
using JLD
include("../src/config.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--scratch"
        help = "Path to scratch directory"
        arg_type = String
        default = "/scratch/gailin/" # local scratch
    "--perplex", "-p"
        help = "Path to PerpleX directory (data files and utilities)"
        arg_type = String
        required = true
    "--perplex_dataset"
    	help = "Which perplex thermodynamic dataset to use? eg hpha02ver.dat or hp02ver.dat"
    	arg_type = String
    	default = "hpha02ver.dat"
end 
parsed_args = parse_args(ARGS, s)
perplex = parsed_args["perplex"]
scratch = parsed_args["scratch"]

#solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)"*
#		"\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar"*
#		"\nDo(HP)\n"
npoints = 20
fluid_endmembers = "abL\nanL\ndiL\nenL\nfaL\nfliq\nfoL\nkspL\nmliq\nqL\nsiL\nq8L\nfa8L\nfo8L\nsil8L\nh2oL\nh2o8L\n"

# From brenhin 
solutions = "Augite(G)\nOpx(JH)\ncAmph(G)\noAmph(DP)\nfeldspar_B\nO(JH)\nSp(JH)\nGrt(JH)\nMica(W)\nBio(TCC)"*
"\nChl(W)\nCtd(W)\nCrd(W)\nSa(WP)\nSt(W)\nIlm(WPH)\nAtg(PN)\nT\nB\nF\nDo(HP)\nScap\nChum\n"


excludes = fluid_endmembers*"ged\nfanth\ngl\nilm\nilm_nol\n"


# "SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2"
# From abers hacker 
vals = transpose([51.0610	49.4020	49.8710	48.9	52.6	45.0	46.1	43.9	44.5	45.5	41.0	40.3	42.5	46.4	47.7	49.6	50.9	50.9	49.1	47.2	47.1	47.6	47.1	48.6	47.0	48.6
0.0000	0.0000	0.0000	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	5.1	5.0	5.1	3.0	3.0	4.5	3.0
15.6551	14.1272	16.1371	3.5	0.0	3.7	3.5	0.8	0.0	1.8	0.0	15.6	17.7	17.5	17.0	16.9	14.3	13.7	15.4	15.3	15.2	14.4	15.4	15.6	15.3	15.1
5.3589	8.5590	5.0190	7.0	4.1	8.7	6.4	9.2	8.2	6.9	9.3	13.1	14.2	11.8	11.2	12.6	12.1	13.8	13.7	10.1	10.5	11.3	13.0	12.4	11.5	12.5
7.0895	11.3563	10.2143	26.4	21.6	40.1	39.0	45.3	47.4	43.6	49.7	7.0	6.7	7.3	7.8	7.4	7.3	7.8	7.6	6.3	6.4	6.6	6.9	7.2	6.8	7.4
19.0322	15.5996	17.3532	13.8	21.7	2.6	5.1	0.8	0.0	2.2	0.0	9.1	9.6	10.7	11.3	10.1	11.7	10.8	11.8	9.1	8.8	10.0	11.7	10.8	10.3	11.1
1.6882	0.9596	1.4090	0.3	0.0	0.0	0.0	0.0	0.0	0.0	0.0	2.1	2.1	2.2	2.4	2.3	2.3	2.3	2.2	2.1	2.6	2.4	2.7	2.5	2.6	2.6
0.0000	0.0000	0.0000	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.2	0.2	0.2	0.2	0.0	0.0	0.0	0.0	0.2	0.2	0.2	0.2	0.2	0.2	0.2
0.0000	0.0000	0.0000	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	7.3	4.5	3.3	2.0	1.2	1.3	0.5	0.0	5.6	5.4	3.1	0.6	0.2	3.0	0.1
0.0000	0.0000	0.0000	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	5.2	2.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0])

### Choose some conditions 
# median depth to 550 isotherm 
isotherm = 39.0;
dpdz = 2900. * 9.8 / 1E5 * 1E3; # km/bar
# sample at half of median depth of crust
depths = [18.355/2, 18.355, 18.355*(3/2)]
# calc t and p 
tts = (depths ./ 39.0).*550 # linear from 0 at surface to 550 at isotherm
pps = depths .* dpdz 
geotherm = 550.0/isotherm/dpdz
P_range = [1, depths[2]*2*dpdz] # to bottom of crust

out = fill(0.0, (3*size(vals,1), 3))
out_modes = Dict()

for d in 1:length(depths)
	for i in 1:size(vals,1)
		# Run perplex
		comp = vals[i,:]
		perplex_configure_geotherm(perplex, scratch, comp, PERPLEX_COMPOSITION_ELTS,
		    P_range, 273.15, geotherm, dataset=parsed_args["perplex_dataset"], solution_phases=solutions,
		    excludes=fluid_endmembers, npoints=npoints)
		seismic = perplex_query_seismic(perplex, scratch)
		selected = findfirst(x -> x > pps[d], seismic["P(bar)"])
		if (seismic["P(bar)"][selected] - pps[d]) > (seismic["P(bar)"][selected-1] - pps[d])
			selected = selected-1
		end
		println("elt $i at $(seismic["T(K)"][selected]-273.15) deg C and $((seismic["P(bar)"][selected])) is 
					vp $((seismic["vp,km/s"][selected])), vs $((seismic["vp,km/s"][selected])/(seismic["vp/vs"][selected])), rho $((seismic["rho,kg/m3"][selected]))")
		#out[i,1] = seismic["T(K)"][selected]-273.15
		#out[i,2] = (seismic["P(bar)"][selected]) / 10000
		out_add = (d-1)*size(vals,1);
		out[out_add+i,1] = (seismic["vp,km/s"][selected])
		out[out_add+i,2] = (seismic["vp,km/s"][selected])/(seismic["vp/vs"][selected])
		out[out_add+i,3] = (seismic["rho,kg/m3"][selected])

		modes = perplex_query_modes(perplex, scratch, include_fluid="y")
		selected = findfirst(x -> x > pps[d], modes["P(bar)"])
		out_modes[(i,depths[d])] = Dict()
		for elt in modes["elements"]
			out_modes[(i,depths[d])][elt] = modes[elt][selected]
		end
		print("$(length(keys(out_modes))) in mode dictionary")
	end 
end

other = readdlm("data/abersHackerRes.csv", ',')
out = hcat(other, out)
println(size(out))
header = ["AH vp" "AH vs" "AH rho" "depth" "P vp" "P vs" "P rho"]
out = vcat(header, out)
println(size(out))
writedlm("data/perplexRes.csv", out, ',')



save("data/perplexModes.jld", "mode_dict", out_modes)




