"""
Alternative to full inversion. 

For a perplex result file f, weight each perplex sample as the sum of its 3d gaussian pdf from each 
unique Crust1.0/geotherm bin sample, then sample from those weighted perplex samples. 
"""

struct SampleModel 
	weights::Array{Float64, 1}
	comps::Array{Float64, 2}
end

"""
	SampleModel(crust, comp)

Arguments: 
crust::Array{Float64, 2} with 
	rows     = unique Crust1.0 sample / geotherm bin combinations 
	columns  = geotherm bin, rho, vp, vpvs 
comp::Array{Float64, 3} with 
	rows     = all perplex samples 
	columns  = columns for inversion 
	z = geotherm 
"""
function SampleModel(crust::Array{Float64, 2}, comp::Array{Float64, 3})

end 