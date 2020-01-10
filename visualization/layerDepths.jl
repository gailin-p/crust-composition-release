## Map crust1 layer depths, find average crust1 layer depths

## For all crust1 squares, depth of bottom crust layer is = to TC1 isotherm.


using StatGeochem
using Statistics
using ProgressMeter
using Plots; gr();

nlong = 180*2
nlat = 90*2

thick = Array{Float64,3}(undef,(4,nlong,nlat))
depth = Array{Float64,3}(undef,(4,nlong,nlat))
isotherm = Array{Float64,2}(undef,(nlong,nlat))

#lats = repeat(-90:90, nlong) # -90 through 90, nlong times
longs = Array{Float64,1}()
lats = Array{Float64,1}()

for lat in -89:90
    for long in -180:179
        push!(longs, long)
        push!(lats, lat)
    end
end

# # Check continent
cont = reshape(map(c -> continents[c], find_geolcont(lats,longs)),(nlong,nlat))
test = cont .== "NA"

# cont = map(c -> continents[c], find_geolcont(lats,longs))
#println(length(cont[cont.=="NA"]))
# lats = lats[cont.!="NA"]
# longs = longs[cont.!="NA"]

# for (i,l) in enumerate(-180:180)
#     start = (i-1)*nlat
#     longs[start+1:start+nlat] = fill(l, nlat)
# end

#longs = vcat(fill.(Array(-180:180), fill(nlat,nlong))) # -180 through 180, each repeated nlat times

for (i, layer) in enumerate([5,6,7,8]) # lower seds, upper, middle, lower
    #if continents[find_geolcont(lat,lon)] != "NA"
    t = reshape(find_crust1_thickness(lats, longs, layer),(nlong,nlat))
    t[test] .= NaN
    thick[i,:,:] = t
    d = reshape(find_crust1_base(lats, longs, layer),(nlong,nlat))
    d[test] .= NaN
    depth[i,:,:] = d
    #isotherm[test] .= NaN
    #println(mean(find_crust1_thickness(lats, longs, layer)))
    #println(mean(find_crust1_base(lats, longs, layer)))
end



#println(nanmean(thick[1,:,:]))
#println(nanmean(thick[2,:,:]))
#println(nanmean(thick[3,:,:]))
#println(nanmean(thick[4,:,:]))
println()
println(nanmean(depth[1,:,:]))
println(nanmean(depth[2,:,:]))
println(nanmean(depth[3,:,:]))
println(nanmean(depth[4,:,:]))

# Look at some histogramz
ds = filter(x -> !isnan(x), depth[2,:,:])
#p = histogram(ds)
#savefig(p, "../output/junk/depth_hist.pdf")


# #println(length(cont[cont.=="NA"]))
# lats = lats[cont.!="NA"]
# longs = longs[cont.!="NA"]

# Display in map
#p = heatmap(thick[3,:,:])
#savefig("../output/thickMap.pdf")
