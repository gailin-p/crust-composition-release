"""
Go directly to perplex for seismic adjustment for dif T, P 
Alternative to parsePerplex.jl + seismic.jl 

NEVER MIND this approach doesn't work because not all phases exist and have properties calculated at every T, P 
"""

"""
Details for running perplex. 
If running threaded, each thread should have its own runner with a different index. 
"""
struct PerplexRunner
    perplex::String
    scratch::String
    dfile::String
    index::Integer 
end 

function PerplexRunner(perplex::String, scratch::String, dfile::String)
	return PerplexRunner(perplex, scratch, dfile, 1)
end 

struct PerplexSeismicError <: Exception 
    msg::String 
end

Base.showerror(io::IO, e::SeismicError) = print(io, "Error in perplex or seismic calc: ", e.msg)


function runSample(pr::PerplexRunner, comp::Array, formationT::Number, formationP::Number, targetT::Array, targetP::Array)
	# Calculate pseudosection 
	# build 
	elementstring = join(uppercase.(elements) .* "\n")
	pseudo_input = "$name\n$dataset\nperplex_option.dat\nn\n2\nn\nn\nn$elementstring\n5\n2\n$tmin\n$tmax\n$pmin\n$pmax\ny\n"
	# Whole-rock composition
    for i = 1:length(composition)
        write(fp,"$(composition[i]) ")
    end
    "n\ny\nsolution_model.dat\n$solutions\n"
	# vertex 

	# werami at formation T, P 


	# Get composition breakdown for formation T and P 

	# Query points in target T and P 
end 