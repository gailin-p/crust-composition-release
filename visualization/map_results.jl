using ArgParse
using DelimitedFiles 
using Plots; gr(); 
using Statistics 
using StatGeochem

include("../src/bin.jl")
include("../src/utilities.jl")

s = ArgParseSettings()
@add_arg_table! s begin
    "--data_file"
        help = "Result data file. For a result from inversion_binned_geotherm"
        arg_type= String
        required=true 
    "--output"
    	help = "Output file name/location. If not given, no file saved."
    	arg_type = String
    	default=""
    "--square_size"
        help = "# latitude /longitude degrees along side of averaged squares"
        arg_type = Int
        default=5
end 

parsed_args = parse_args(ARGS, s)

(dat, header) = readdlm(parsed_args["data_file"], ',', header=true)
header = header[:] # make 1d 

si_i = findfirst(isequal("SiO2"), header) # what we'll plot  
lat_i = findfirst(isequal("sample_lat"), header)
long_i = findfirst(isequal("sample_long"), header)

# Take area average 
lats, longs, si = areaAverage(Float64.(dat[:,lat_i]), 
    Float64.(dat[:,long_i]), Float64.(dat[:,si_i]), size=parsed_args["square_size"])
globe_s = globe(lats, longs, si)
println(nanmean(globe_s))
#globe_s[1,1] = 85; globe_s[1,2] = 40 # put on same scale 

#p = plot(longs, lats, seriestype=:scatter, markerstrokewidth=0, marker_z=si, colorbar=true)
p1 = heatmap(globe_s, grid=false, axis=false, title="Inverted composition", 
    colorbar_title="% SiO2", c=:haline);

display(p1)

## Compare to area averaged compositions 
(dat_est, header_est) = readdlm("data/remote/base/bsr_ignmajors_1.csv", ',', header=true)
header_est = header_est[:]
lat_i = findfirst(isequal("Latitude"), header_est)
long_i = findfirst(isequal("Longitude"), header_est)
si_i = findfirst(isequal("SiO2"), header_est)

s_lats, s_longs, s_si = areaAverage(Float64.(dat_est[:,lat_i]), 
    Float64.(dat_est[:,long_i]), Float64.(dat_est[:,si_i]), size=parsed_args["square_size"])

surface_globe = globe(s_lats, s_longs, s_si)
p_surface = heatmap(surface_globe, grid=false, axis=false, title="Surface composition", 
    colorbar_title="% SiO2", c=:haline);

m = Dict{Tuple, Tuple}() # map from (lat,long) to (estimated,exposed)
for i in 1:length(lats)
	m[(lats[i],longs[i])] = (si[i], NaN)
end 
for i in 1:length(s_lats)
	v = get!(m, (s_lats[i], s_longs[i]), (NaN,NaN))
	m[(s_lats[i], s_longs[i])] = (v[1], s_si[i])
end

estimates = [v[1] for v in values(m)]
surface = [v[2] for v in values(m)]
p2 = scatter(surface, estimates, markerstrokewidth=0, markeralpha=.3, ylim=(40,80), xlim=(40,80),
	legend=false, xlabel="Area-averaged exposed composition (% SiO2)", 
	ylabel="Estimated composition (% SiO2)");

x, y, yerr = bin(surface, estimates, 40, 80, 10)
plot!(p2, x, y, yerror=yerr);

p = plot(p1, p_surface, p2, layout=(3,1), size=(400,700))
display(p)

if length(parsed_args["output"])>0
	savefig(p, parsed_args["output"])
else 
	sleep(60*4) # give user a chance to look at plot 
end 













