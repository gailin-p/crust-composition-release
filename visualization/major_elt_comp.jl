"""
Given the final run output (csv saved by inversion_binned_geotherm), 
give the major elt composition for all layers. 
"""

using ArgParse

include("../src/config.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/bin_geotherm_si_weighted"
    "--model", "-m"
        help = "Type of model to use. Allowed: range. Inversion model doesn't save indices currently."
        arg_type = String
        range_tester = x -> (x in ["range"])
        default = "range"
end
parsed_args = parse_args(ARGS, s)

function run(parsed_args)
	# Read inversion file 
	filename = "../data/"*parsed_args["data_prefix"]*"/inverted_samples_$(parsed_args[model]).csv"
	indices, header = readdlm(filename, ',', header=true)

	# Read composition file. Compositions are same for all geotherms 
	filename = "data/"*parsed_args["data_prefix"]*"/bsr_ignmajors_1.csv" 
	compositions = readdlm(filename, ',')

	# For each major elt...
	for elt_name in COMPOSITION_ELEMENTS
		elt_i = findfirst(isequal(elt_name), PERPLEX_ELEMENTS)
	# For each layer, take average of samples in each row, then average of those averages. 


	# todo: Propogate error sensibly 


end 

run(parsed_args)