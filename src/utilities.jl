using Statistics
using DelimitedFiles 
using StatGeochem 
using HDF5

"""
Keep a running mean
m mean so far 
n number of samples so far 
"""
mutable struct RunningMean
	m::Float64
	n::Int
end 

function RunningMean()
	return RunningMean(0., 0)
end 

"""
Add a new number to a running mean 
"""
function mean!(running::RunningMean, new::Number)
	running.m = (running.n*running.m + new)/(running.n + 1)
	running.n += 1 
end 

"""
Normalize each row of input matrix (in place)
"""
function normalizeComp!(a::AbstractArray{Float64, 2})
	sums = sum(a,dims=2)
	for row in 1:size(a,1)
		a[row,:] .= (a[row,:] ./ sums[row]) .* 100
	end 
end 

function inverseMean(a::AbstractArray{Float64})
	return 1/(nanmean(1 ./ a))
end 

"""
Write options to option file
"""
function writeOptions(filename, options)
	out = fill("", (length(keys(options)), 2))
	for (k, key) in enumerate(keys(options))
		out[k, 1] = key 
		out[k, 2] = string(options[key])
	end 
	writedlm(filename, out, ",")
end 

"""
For values at each location in latitide/longitude, find average value per lat/long grid square.
Return lat, long, ave value for each square, with val=NaN for squares with no values. 
size combines that size^2 lat/long squares into each returned square. 
"""
function areaAverage(latitude::Array{Float64,1}, longitude::Array{Float64,1}, vals::Array{Float64,1}; size::Int=1)
	good = .!(isnan.(latitude) .| isnan.(longitude) .| isnan.(vals))

	m = Dict{Tuple, Array}() # map from lat/long to list of values
	# Want full globe for visualization, even where no values 
	# for lat in range(floor(minimum(latitude[good])/size), stop=floor(maximum(latitude[good])/size))
	# 	for long in range(floor(minimum(longitude[good])/size), stop=floor(maximum(longitude[good])/size))
	# 		m[(floor(lat), floor(long))] = [] 
	# 	end
	# end

	# now put in those values 
	latitude = latitude[good] ./ size 
	longitude = longitude[good] ./ size
	vals = vals[good]
	for i in 1:length(latitude)
		a = get!(m, (floor(latitude[i]), floor(longitude[i])), [])
		append!(a, vals[i])
	end 
	for k in keys(m)
		m[k] = [nanmean(Array{Float64}(m[k]))]
	end 

	lats = [k[1] for k in keys(m)]
	longs = [k[2] for k in keys(m)]
	val = [mean(v) for v in values(m)]

	#good = .!(isnan.(lats) .| isnan.(longs) .| isnan.(val))

	#return lats[good].*size, longs[good].*size, val[good]
	return lats.*size, longs.*size, val
end 

"""
For plotting latitude and longitude data.
Assumes lat/long pairs are unique (run after areaAverage)
Returns grid with NaN at any missing values 
"""
function globe(lats, longs, vals)
	slats = sort(unique(lats))
	slongs = sort(unique(longs))
    globe = fill(NaN, (length(slats), length(slongs)))
    for i in 1:length(lats)
        y = searchsortedfirst(slats, lats[i])
        x = searchsortedfirst(slongs, longs[i])
        #println("$x, $y")
        globe[y,x] = vals[i]
    end 
    return globe 
end 

function plotglobe(k::String, dat, header)
    return globe( areaAverage( 
    	dat[:, findfirst(isequal("sample_lat"), header)], 
    	dat[:, findfirst(isequal("sample_long"), header)], 
    	dat[:, findfirst(isequal(k), header)])...)
end

"""
Take area averages of composition samples and corresponding Perple_X samples. 
"""
function areaAverage(inputDataDir::String, outputDataDir::String)
	# use file from first geotherm bin to collect indices by lat/long 
	file = "data/" * inputDataDir * "/bsr_ignmajors_1.csv"
	dat = readdlm(file, ',')
	mkdir("data/$outputDataDir")

	# Create mapping from unique lat/longs to sample indices 
	lat_i = findfirst(isequal("Latitude"), PERPLEX_ELEMENTS)
	long_i = findfirst(isequal("Longitude"), PERPLEX_ELEMENTS)
	dat[:,lat_i] = floor.(dat[:,lat_i])
	dat[:,long_i] = floor.(dat[:,long_i])
	m = Dict{Tuple, Array}()
	for j in 1:size(dat,1)
		a = get!(m, (dat[j, lat_i], dat[j, long_i]), [])
		append!(a, j)
	end 

	# Average data by lat/long bin for every (geotherm bin)/layer combo 
	fileNames = filter(x->contains(x,"perplex_out_"), readdir("data/$inputDataDir"))
	nBins = length(fileNames)
	for file_i in 1:nBins 
		ign = readdlm("data/$inputDataDir/bsr_ignmajors_$file_i.csv", ',')
		perplex = h5read("data/$inputDataDir/perplex_out_$file_i.h5", "results")
		new_perplex = fill(NaN, (4,3,length(m)))
		new_ign = fill(NaN, (length(m), size(dat,2)))
		for (new_i, k) in enumerate(keys(m))
			new_ign[new_i, :] = mean(dat[m[k],:], dims=1)
			new_ign[new_i, 1] = new_i
			for j in 1:size(perplex, 1)
				for l in 1:size(perplex, 2)
					new_perplex[j, l, new_i] = inverseMean(perplex[j, l, m[k]])
				end 
			end 
			new_perplex[1, :, new_i] .= new_i # index is not average index 
		end 
		writedlm("data/$outputDataDir/bsr_ignmajors_$file_i.csv", new_ign, ',')
		h5write("data/$outputDataDir/perplex_out_$file_i.h5", "results", new_perplex)
	end 

	# Write options file 
	writeOptions("data/$outputDataDir/areaAverage_options.csv", Dict([("inputDir", inputDataDir)]))
end 














