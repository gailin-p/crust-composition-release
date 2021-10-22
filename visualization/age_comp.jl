"""
Create main figure (composition of upper, middle, and lower crust over time)

"""

using ArgParse
using DelimitedFiles
using Plots; gr();
using Statistics

include("../src/bin.jl")
include("../src/config.jl")
include("../src/crustDistribution.jl")
include("../src/inversionModel.jl")
include("../src/invertData.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--data_path", "-p"
        help = "Path to model output"
        arg_type= String
        required=true
end

parsed_args = parse_args(ARGS, s)


build_args = readdlm("$(parsed_args["data_path"])/inversion_options.csv", ',', header=false)
N = build_args[findfirst(isequal("num_invert"), build_args[:,1]),2]
M = build_args[findfirst(isequal("num_runs"), build_args[:,1]),2]

# By age!
function find_age_aves(files, N, M)
	age_model = EarthChemAge(10, 3) # Default in inversion_binned_geotherm
	ages = ageBins(age_model)


    age_median = Array{Float64, 2}(undef, 3, length(ages)-1) # (layer, age)
	age_low = Array{Float64, 2}(undef, 3, length(ages)-1) # (layer, age)
	age_high = Array{Float64, 2}(undef, 3, length(ages)-1) # (layer, age) # 95 percentile


	for (l, file) in enumerate(files)
		layer = LAYER_NAMES[l]

		(result, header) = readdlm(file, ',', header=true)
		header = header[:] # make 1d

		si_index = findfirst(isequal("SiO2"), header) # result
		age_index = findfirst(isequal("sample_age"), header) # Age

		for (i, age) in enumerate(ages[1:end-1])
			mc_avgs = fill(NaN, M)
			for m in 1:M
				# only this run
				istart = (m-1)*N + 1
				iend = m*N

				dat = result[istart:iend, :]

                good = ((.~isnan.(dat[:,si_index])) .& (.~isnan.(dat[:,age_index])))
				test = (dat[:,age_index] .>= age) .& (dat[:,age_index] .< ages[i+1]) .& good
                mc_avgs[m] = nanmean(convert(Array{Float64,1}, dat[test,si_index]))
				#println("$(sum(test)) points of age $age")
            end

            check = .! isnan.(mc_avgs)
			try
            	age_median[l,i] = median(mc_avgs[check])
            	age_low[l,i] = percentile(mc_avgs[check], 25)
            	age_high[l,i] = percentile(mc_avgs[check], 75)
			catch ArgumentError # if empty array
				age_median[l,i] = NaN
            	age_low[l,i] = NaN
            	age_high[l,i] = NaN
			end
		end
	end
	return ages, age_median, age_low, age_high
end

files = ["$(parsed_args["data_path"])/results-$layer.csv" for layer in LAYER_NAMES]
ages, age_results, age_low, age_high = find_age_aves(files, N, M)


p1 = plot(legend=false, ylabel="% SiO2", ylims=(52, 70), title="Upper", titlefontsize=11)
p2 = plot(legend=:bottomright, legendfontsize=7, fg_legend = :transparent, ylabel="% SiO2", ylims=(52, 70), title="Middle", titlefontsize=11)
p3 = plot(legend=false, xlabel="Age", ylabel="% SiO2", ylims=(52, 70), title="Lower", titlefontsize=11)
ps = [p1, p2, p3]

p = plot(size=(550,400), legend=:bottomright, fg_legend = :transparent, framestyle=:box, xlabel="Age", ylabel="SiO2 (weight %)");

for i in 1:3 # layers
	color = [:blue, :orange, :green][i]
    plot!(p, ages[1:end-1].+i*10, age_results[i,:],
        ribbon=(age_results[i,:].-age_low[i,:], age_high[i,:].-age_results[i,:]),
        label=LAYER_NAMES[i], markerstrokecolor=color, fillalpha=.35,
        marker=true, markersize=3,
		linecolor=color, markercolor=color)
end



# Plot exposed. Always use base so have same comparison...
# tmin = 0
# tmax = 4000
# nbins = 4
# samplePath = "data/remote/base_nobin/bsr_ignmajors_1.csv"
# ign, h = readdlm(samplePath, ',', header=true)
# age_i = findfirst(isequal("Age"),PERPLEX_ELEMENTS)
# si_i = findfirst(isequal("SiO2"),PERPLEX_ELEMENTS)
# original = matread(IGN_FILE)
# age_centers, elt_means, elt_errors = bin(ign[:,age_i], ign[:,si_i],
#         tmin, tmax, length(ign[si_i])/length(original["SiO2"]), nbins)
#
# plot!(ps[1], age_centers, elt_means, color=:pink,
# 	yerror=elt_errors, label="exposed", markerstrokecolor=:auto);
# plot!(p, age_centers, elt_means, color=:pink,
# 	yerror=elt_errors, label="exposed", markerstrokecolor=:auto);


outputPath = "$(parsed_args["data_path"])/output"
mkpath(outputPath) # make if does not exist
savefig(p, outputPath*"/composition-ages.pdf");
