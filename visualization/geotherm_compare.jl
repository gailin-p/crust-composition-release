using ArgParse
using DelimitedFiles 
using Plots; gr(); 
using Statistics 
using StatGeochem

include("../src/bin.jl")
include("../src/utilities.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_file"
        help = "Result data file. For a result from inversion_binned_geotherm"
        arg_type= String
        required=true 
    "--output"
    	help = "Output file name/location. If not given, no file saved."
    	arg_type = String
    	default=""
end 

parsed_args = parse_args(ARGS, s)

(dat, header) = readdlm(parsed_args["data_file"], ',', header=true)
header = header[:] # make 1d 

println(header)

si_i = findfirst(isequal("SiO2"), header) # what we'll plot  
bin_i = findfirst(isequal("bin"), header)
geo_i = findfirst(isequal("geotherm"), header)

lat_i = findfirst(isequal("sample_lat"), header)
long_i = findfirst(isequal("sample_long"), header)

# Take area average 
lats, longs, si = areaAverage(dat[:,lat_i], dat[:,long_i], dat[:,si_i])
globe_si = globe(lats, longs, si)

# Geotherm average 
lats, longs, gbin = areaAverage(dat[:,lat_i], dat[:,long_i], dat[:,bin_i])
globe_bin = globe(lats, longs, gbin)

# Geotherm of samples 
lats, longs, geotherms = areaAverage(dat[:,lat_i], dat[:,long_i], dat[:,bin_i])
globe_geo = globe(lats, longs, geotherms)

# Average geotherm of 


#p = plot(longs, lats, seriestype=:scatter, markerstrokewidth=0, marker_z=si, colorbar=true)
p1 = heatmap(globe_si, axis=false, title="SiO2");
p2 = heatmap(globe_bin, axis=false, title="Bin used");
p3 = heatmap(globe_geo, axis=false, title="Geotherm of sample");

p = plot(p1, p2, p3)
display(p)

sleep(60)



