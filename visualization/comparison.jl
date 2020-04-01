"""
Compare our results to results of other workers 
"""

mutable struct Result
	name::String
	composition::Array{Float64, 1} # Lower, middle, upper
	#depth::Array{Float64, 1} # Depths for estimates 
	notes::String
end 

results = [
	Result("Rudnick and Fountain", [52.3, 60.6, 66.0], "upper crust from Taylor and McLennan 1985"),
	Result("Rudnick and Gao", [53.4, 63.5, 66.6], ""),
	Result("Hacker, Kelemen, Behn (most felsic)", [-1, 68, 62], 
		"Most mafic felsic model of the three-layer models in this paper.")
	Result("Hacker, Kelemen, Behn (most mafic)", [-1, 53, 49],
		"Most mafic model of the three-layer models in this paper.")
	
]