## ---
    using StatGeochem, Plots
    using NetCDF

## --- Shan 2016

    ncinfo("data/US.2016.nc")

    # Dims: longitude latitude depth. NaN-value = 9999.0
    vp_us = ncread("data/US.2016.nc", "vp") .|> x -> x==9999.0 ? NaN : x
    vs_us = ncread("data/US.2016.nc", "vsv") .|> x -> x==9999.0 ? NaN : x
    rho_us = ncread("data/US.2016.nc", "rho") .|> x -> x==9999.0 ? NaN : x

    depth_us = ncread("data/US.2016.nc", "depth")
    latitude_us = ncread("data/US.2016.nc", "latitude")
    # values range, bizarely, from 235 to 295. This puts it in the right range but idk if correct. 
    longitude_us = ncread("data/US.2016.nc", "longitude") .- 360 

## --- Plot a random selection of Vp-depth paths (Shan)

    h = plot(xlabel="Vp", ylabel="Depth")
    for i=1:500
        lat = rand(1:length(latitude_us))
        lon = rand(1:length(longitude_us))
        plot!(h, vp_us[lon, lat, :] |> vec, -depth_us, label="", color=:darkblue, alpha=0.2)
    end
    display(h)

## --- Plot a random selection of Vs-depth paths (Shan)

    h = plot(xlabel="Vs", ylabel="Depth")
    for i=1:500
        lat = rand(1:length(latitude_us))
        lon = rand(1:length(longitude_us))
        plot!(h, vs_us[lon, lat, :] |> vec, -depth_us, label="", color=:darkblue, alpha=0.2)
    end
    display(h)

## --- Plot a random selection of rho-depth paths (Shan)

    h = plot(xlabel="Density", ylabel="Depth")
    for i=1:500
        lat = rand(1:length(latitude_us))
        lon = rand(1:length(longitude_us))
        plot!(h, rho_us[lon, lat, :] |> vec, -depth_us, label="", color=:darkblue, alpha=0.2)
    end
    display(h)

## --- Screen out mantle values (Shan)

    ismantle = (rho_us .> 3.2) .| (vp_us .> 7.2) .| (repeat(reshape(depth_us, (1,1,length(depth_us))), outer=(size(rho_us)[1:2]...,1)) .> 70)
    rho_crust_us = rho_us .+ ismantle.*NaN # preserves dimensions 
    vp_crust_us = vp_us .+ ismantle.*NaN
    vs_crust_us = vs_us .+ ismantle.*NaN


    vp_crust_ave_us = nanmean(vp_crust_us, dim=(1,2))
    vs_crust_ave_us = nanmean(vs_crust_us, dim=(1,2))
    rho_crust_ave_us = nanmean(rho_crust_us, dim=(1,2))


## --- Wang, 2020, Southern California

    # ncinfo("SoCal.ANAT-Vs+RA.Wang.2020.nc")

    # # Dims: longitude latitude depth. NaN-value = -1000.0
    # vp_sc = ncread("SoCal.ANAT-Vs+RA.Wang.2020.nc", "voigtVp") .|> x -> x==-1000.0 ? NaN : x
    # vs_sc = ncread("SoCal.ANAT-Vs+RA.Wang.2020.nc", "voigtVs") .|> x -> x==-1000.0 ? NaN : x
    # rho_sc = ncread("SoCal.ANAT-Vs+RA.Wang.2020.nc", "rho") .|> x -> x==-1000.0 ? NaN : x

    # depth_sc = ncread("SoCal.ANAT-Vs+RA.Wang.2020.nc", "depth")
    # latitude_sc = ncread("SoCal.ANAT-Vs+RA.Wang.2020.nc", "latitude")
    # longitude_sc = ncread("SoCal.ANAT-Vs+RA.Wang.2020.nc", "longitude")

## ---  Plot a random selection of Vp-depth paths (Wang)

#     h = plot(xlabel="Vp", ylabel="Depth")
#     for i=1:500
#         lat = rand(1:length(latitude_sc))
#         lon = rand(1:length(longitude_sc))
#         plot!(h, vp_sc[lon, lat, :] |> vec, -depth_sc, label="", color=:darkblue, alpha=0.2)
#     end
#     display(h)

# ## ---  Plot a random selection of Vs-depth paths (Wang)

#     h = plot(xlabel="Vs", ylabel="Depth")
#     for i=1:500
#         lat = rand(1:length(latitude_sc))
#         lon = rand(1:length(longitude_sc))
#         plot!(h, vs_sc[lon, lat, :] |> vec, -depth_sc, label="", color=:darkblue, alpha=0.2)
#     end
#     display(h)

# ## --- Plot a random selection of rho-depth paths (Wang)

#     h = plot(xlabel="Density", ylabel="Depth")
#     for i=1:200
#         lat = rand(1:length(latitude_sc))
#         lon = rand(1:length(longitude_sc))
#         plot!(h, rho_sc[lon, lat, :] |> vec, -depth_sc, label="", color=:darkblue, alpha=0.2)
#     end
#     display(h)

# ## ---  Screen out mantle values (Wang)

#     ismantle = (rho_sc .> 3.1) .| (vp_sc .> 7.2) .| (repeat(reshape(depth_sc, (1,1,length(depth_sc))), outer=(size(rho_sc)[1:2]...,1)) .> 40)
#     rho_crust_sc = rho_sc .+ ismantle.*NaN
#     vp_crust_sc = vp_sc .+ ismantle.*NaN
#     vs_crust_sc = vs_sc .+ ismantle.*NaN

#     vp_crust_ave_sc = nanmean(vp_crust_sc, dim=(1,2))
#     vs_crust_ave_sc = nanmean(vs_crust_sc, dim=(1,2))
#     rho_crust_ave_sc = nanmean(rho_crust_sc, dim=(1,2))


## --- Crust 1.0

    using MAT
    igncn1 = matread("igncn1.mat")


## --- Distribute over width of crustal layer

    Depth = [-igncn1["Upper_Crust"].*rand.();
        -igncn1["Upper_Crust"]-igncn1["Middle_Crust"].*rand.();
        -igncn1["Upper_Crust"]-igncn1["Middle_Crust"]-igncn1["Lower_Crust"].*rand.() ]

    Vp = [igncn1["Upper_Vp"]; igncn1["Middle_Vp"]; igncn1["Lower_Vp"]] .*(1 .+ 0.025.*randn.())
    Vs = [igncn1["Upper_Vs"]; igncn1["Middle_Vs"]; igncn1["Lower_Vs"]] .*(1 .+ 0.025.*randn.())
    Rho = [igncn1["Upper_Rho"]; igncn1["Middle_Rho"]; igncn1["Lower_Rho"]] .*(1 .+ 0.025.*randn.())

## ---

    (c,m,e) = binmeans(Depth |> vec, Vp |> vec, -80, 0, 20)
    h = plot(m,c,xerror=e, xlabel="Vp (km/s)", ylabel="depth (km)", label="Crust 1.0")
    plot!(h,vp_crust_us,-depth_us, label="Wang 2016 (US)")
    # plot!(h,vp_crust_ave_sc,-depth_sc, label="Shan 2020 (So-Cal)")
    plot!(h,xlims=(4.5,7.2),ylims=(-70,0), fg_color_legend=:white, framestyle=:box, legend=:bottomleft)
    savefig(h,"Vp_comparison.pdf")
    display(h)

## ---

    (c,m,e) = binmeans(Depth |> vec, Vs |> vec, -80, 0, 20)
    h = plot(m,c,xerror=e, xlabel="Vs (km/s)", ylabel="depth (km)", label="Crust 1.0")
    plot!(h,vs_crust_us,-depth_us, label="Wang 2016 (US)")
    # plot!(h,vs_crust_ave_sc,-depth_sc, label="Shan 2020 (So-Cal)")
    plot!(h,xlims=(2,4.2),ylims=(-70,0), fg_color_legend=:white, framestyle=:box, legend=:bottomleft)
    savefig(h,"Vs_comparison.pdf")
    display(h)

## ---

    (c,m,e) = binmeans(Depth |> vec, Rho |> vec, -80, 0, 20)
    h = plot(m,c,xerror=e, xlabel="rho (kg/m3)", ylabel="depth (km)", label="Crust 1.0")
    plot!(h,rho_crust_us,-depth_us, label="Wang 2016 (US)")
    # plot!(h,rho_crust_ave_sc,-depth_sc, label="Shan 2020 (So-Cal)")
    # plot!(h,xlims=(2,3.1),ylims=(-70,0), fg_color_legend=:white, framestyle=:box, legend=:bottomleft)
    savefig(h,"rho_comparison.pdf")
    display(h)

## ---
