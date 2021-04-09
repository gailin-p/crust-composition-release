using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "--data_prefix", "-d"
        help = "Folder for data files"
        arg_type= String
        default="remote/base"
    "--num_invert", "-n" 
    	help = "How many crust samples to invert?"
    	arg_type = Int
    	default = nOriginal()
end 


