"""
Code for transforming perplex seismic output to different T and P using properties at calculated T&P 
"""


struct SeismicTransform
    datafile::String
    dat::Dict{String,Dict{String,Any}}
    Tref::Float64
    Pref::Float64
end 

struct SeismicError <: Exception 
    msg::String 
end

Base.showerror(io::IO, e::SeismicError) = print(io, "Error in seismic calc: ", e.msg)


function SeismicTransform(datafile::String)
    return SeismicTransform(datafile, parsePerplexData(datafile), 298.15, 1.0)
end

# Default. TODO include this file in data dir 
function SeismicTransform()
    return SeismicTransform("/users/gailin/resources/perplex-stable/hpha11ver.dat")
end


function get_seismic(T, P, properties, endmembers, st)
    for e in keys(endmembers)
        if !(e in keys(st.dat))
            println("missing $e")
        end
    end

    # Get shear and bulk moduli at new T, P
    # don't include liquid water or co2 in seismic calc 
    delete!(endmembers, "H2O")
    delete!(endmembers, "CO2")

    e_shear = Dict()
    for e in keys(endmembers)
        e_shear[e] = adiabatic_shear(e, T, P, st)
    end 
    
    # Take VRH average of moduli 
    mu = vrh(endmembers, e_shear)
    
    # Get new bulk modulus: linear approx based on T, P derivatives 
    ks = properties["Ks(bar)"] + 
        properties["Ks_P"]*(P-properties["P(bar)"]) + 
        properties["Ks_T(bar/K)"]*(T-properties["T(K)"])
    
    # Get new density 
    rho = density_adjust(T, P, properties["T(K)"], properties["P(bar)"], 
        properties["Alpha(1/K)"], properties["Beta(1/bar)"], properties["Density(kg/m3)"])
    
    # Convert bar -> Pa 
    mu_p = mu*100000  
    k_p = ks*100000
    
    # Get new vp, vs 
    vp = sqrt((k_p + 4*mu_p/3)/rho)/1000 # m/s to km/s
    vs = sqrt(mu_p/rho)/1000 # m/s to km/s
    
    return rho, vp, vs 
end

function adiabatic_shear(endmember, Tnew, Pnew, st)
    try 
        c = st.dat[endmember]
        mu = c["m0"] + c["m1"]*(Pnew - st.Pref) + c["m2"]*(Tnew - st.Tref)
        return mu
    catch e 
        throw(SeismicError("Did not find endmember $endmember"))
    end 
end 

# given vp, vs in km/s, rho in kg/m^3, return mu and ks in bar 
function ks_mu(rho, vp, vs)
    mu = rho*(1000*vs)^2
    ks = rho*(1000*vp)^2 - (4/3)*mu 
    return ks/100000, mu/100000 # convert from Pa to Bar 
end 


"""
    vrh(endmembers, shears)
"""
function vrh(endmembers::Dict, shears::Dict)
    if keys(endmembers) != keys(shears)
        throw(SeismicError("Keys of endmembers and keys of shear velocities do not match"))
    end 
    total = sum(values(endmembers))
    voigt = sum([shears[e]*(endmembers[e]/total) for e in keys(endmembers)])
    reuss = 1/sum([(1/shears[e])*(endmembers[e]/total) for e in keys(endmembers)])
    #println("$voigt, $reuss")
    return (voigt + reuss)/2
end

function density_adjust(Tnew, Pnew, Told, Pold, alpha, beta, density)
    # Adjust assuming unit mass, so V = 1/density, density = 1/volume
    #Vnew = (1/density)+alpha*(Tnew-Told)+beta*(Pnew-Pold)
    V = 1/density
    #println("temp adjustment $(alpha*(Tnew-Told)), pressure adjustment $(beta*(Pnew-Pold))")
    Vnew = V + V*(alpha*(Tnew-Told)) - V*(beta*(Pnew-Pold))
    #println("old $density, new $(1/Vnew)")
    return 1/Vnew
end




