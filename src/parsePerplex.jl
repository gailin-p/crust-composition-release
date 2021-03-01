"""
Utilities for parsing perplex output and data files. 
Intended use: 
 - Get endmember fractions at some rock formation T&P from Perple_X 
 - Get endmember properties (m0, m1, m2, k0, k1, k2) from perplex data file
 - Calculate rock properties (Vs, Rho, Vp) at a new T, P given those properties 
"""

struct ParsePerplexError <: Exception 
    msg::String 
end

# Map from perplex endmembers in solution_model.dat to fundamental endmembers in hphaetc.dat 
pmap = Dict([("om", ["jd", "di"]), # om = 1/2 jd + 1/2 di 
        ("opx", ["en", "fs"]), 
        ("cfm", ["di", "hed"]),
        ("jac", ["jd", "acm"]),
        ("gl_dqf", ["fgl"]),
        ("fanth_dq", ["fanth"]),
        ("ammo1", ["anth", "fanth"]), # actually a 3/7, 4/7 split but too small an effect to matter 
        ("ammo2", ["anth", "fanth"]), # 2/7, 5/7 
        ("ogl_dqf", ["gl"]),
        ("ged_dqf", ["ged"]),
        ("mpa", ["parg", "anth"]), # actually mpa      = 1 parg - 1 tr +1 anth 
        ("omrb_dqf", ["gl", "acm"]), # actually omrb_dqf = 1 gl -2 jd +2 acm
        ("tbit", ["phl", "ru"]), # actually tbit = 1 phl - 1 br + 1 ru 
        ("obi", ["phl", "ann"]) # obi = 2/3 phl + 1/3 ann
        ]); 

Base.showerror(io::IO, e::ParsePerplexError) = print(io, "Error when parsing perplex: ", e.msg)  

function parsePerplexData(dataFile::String)
	if !isfile(dataFile)
		throw(ParsePerplexError("cannot read file $dataFile"))
	end
	data = open(f->read(f, String), dataFile)
	# split on end surrounded by line break and/or space 
	# First is header, last is surrounding end (no elt there)
	elts = split(data, r"[\r\n]+\s*end\s*[\r\n]+")[2:end-1] 
	library = Dict()
	for e in elts
		#println(e)
		try 
			d = parsePerplexItem(e)
			library[d["name"]] = d
		catch err 
			println(err)
            throw(err)
		end
	end 
	return library
end

"""
	parsePerplexItem(item)
Parse one perplex data file entry. return dictionary of features (EoS, H, G0, etc)
"""
function parsePerplexItem(item)
	d = Dict()
	for m in eachmatch(r"[a-zA-Z0-9]+\s*=\s*[a-zA-Z0-9.+-]+", item)
	    parts = strip.(split(m.match, '='))
	    num = replace(lowercase(parts[2]), "d"=>"e") # some entries use d for scientific notation
	    d[parts[1]] = parse(Float64, num)
	end
	name = split(item)[1]
	d["name"] = name 
	return d 
end 


"""
For a perplex point output, return dictionary from endmember phases to the vol% of that phase 
"""
function parse_perplex_point(point)
    titles = [m.match for m in eachmatch(r"[\r\n]+[A-Za-z0-9\s)(/]+:[\r\n]+", point)]
    out = split(point, r"[\r\n]+[A-Za-z0-9\s)(/]+:[\r\n]+")[2:end] # first is space above first title
    toreturn = Dict()
    
    # Parse phase compositions
    phase_comp_i = findfirst(a->occursin("phase compositions", a), lowercase.(titles))
    rows = split(out[phase_comp_i], r"[\r\n]+")
    header = split(rows[1], r"\s{2,}")
    names = [strip(split(r, r"\s{2,}")[1]) for r in rows[2:end]]
    volp = zeros(length(names))
    target = findfirst(isequal("vol %"), header)
    for (i,r) in enumerate(rows[2:end])
        volp[i] = parse(Float64, split(r, r"\s{2,}")[target])
    end
    
    # Parse molar percents 
    molar_i = findfirst(a->occursin("phase speciation", a), lowercase.(titles))
    phases = split(out[molar_i], r"[\r\n]+")

    # Check if any of the phases are just overflow of the last phase: if no title set apart by >2 spaces
    solutions = [strip(split(r, r"\s{2,}")[1]) for r in phases]
    while "" in solutions
        i = findfirst(a->isequal("",a), solutions)
        phases[i-1] = phases[i-1]*", "*strip(phases[i])
        deleteat!(phases, i)
        solutions = [strip(split(r, r"\s{2,}")[1]) for r in phases]
    end 

    breakdown = [strip(split(r, r"\s{2,}")[2]) for r in phases]
    
    # For phases that are not solutions, add directly to toreturn 
    endmembers = [n for n in names if !(n in solutions)]
    if length(endmembers) != length(unique(endmembers))
        throw(ParsePerplexError("Warning! Endmember seems to be included twice in phase composition list. Unable to parse \r\n $point"))
    end
    for e in endmembers
        voli = findfirst(isequal(e), names)
        toreturn[e] = volp[voli]
    end
    
    # solutions will be both in phase composition list and phase speciation list. 
    # For those, iterate through their endmembers, multiplying their volume % by endmember % 
    for (namei, name) in enumerate(names)
        if length(solutions) == 0
            break
        elseif solutions[1] != name
            throw(ParsePerplexError("Warning! Solution order does not match names from phase composition list. Unable to parse \r\n $point."))
        else
            b = popfirst!(breakdown)
            s = popfirst!(solutions)
            for item in split(b, ",")
                n = strip(split(item, ":")[1])
                p = parse(Float64, split(item, ":")[2])
                a = get(toreturn, n, 0.0)
                toreturn[n] = a + p*volp[namei]
            end

        end
    end

    # Use common names 
    toreturn_final = Dict() 
    for name in keys(toreturn) 
        breakdown = get(pmap, name, [name]) 
        for n in breakdown # new names 
            existing = get(toreturn_final, n, 0.0)
            toreturn_final[n] = existing + toreturn[name]/length(breakdown)
        end
    end 

    # No negative endmembers 
    for name in keys(toreturn_final)
        if toreturn_final[name] < 0
            toreturn_final[name] = 0 # some may be < 0 b/c weird perplex endmember math 
        end
    end

    if abs(sum(values(toreturn_final))-100) > 1
        throw(ParsePerplexError("Endmembers do not sum to 1"))
    end 

    return toreturn_final
end 

"""
For point output of perplex, return T, P, density, alpha, beta, Ks, Mu, vp, vs 
"""
function get_system_props(point)
    titles = [m.match for m in eachmatch(r"[\r\n]+[A-Za-z0-9\s)(/]+:[\r\n]+", point)]
    out = split(point, r"[\r\n]+[A-Za-z0-9\s)(/]+:[\r\n]+")[2:end] # first is space above first title
    toreturn = Dict()
    
    # Get alpha, beta, density 
    properties_i = findfirst(a->occursin("properties and density", a), lowercase.(titles))
    if isnothing(properties_i)
        throw(ParsePerplexError("Cannot find properties in perplex output \r\n $point"))
    end
    rows = split(out[properties_i], r"[\r\n]+")
    header = split(rows[1], r"\s{2,}")
    rows = rows[2:end]
    cols = [strip(split(r, r"\s{2,}")[1]) for r in rows]
    if !("System" in cols)
        throw(ParsePerplexError("System not found in molar properties, cannot parse point \r\n $point"))
    end
    systemi = findfirst(isequal("System"), cols)
    targets = ["Alpha(1/K)", "Beta(1/bar)", "Density(kg/m3)"]
    for t in targets
        if !(t in header)
            throw(ParsePerplexError("$t not found in molar properties, cannot parse point \r\n $point"))
        end
        ti = findfirst(isequal(t), header)
        toreturn[t] = parse(Float64, split(rows[systemi], r"\s{2,}")[ti])
    end
    
    if toreturn["Alpha(1/K)"]> 1/100
        throw(ParsePerplexError("Alpha = $(toreturn["Alpha(1/K)"]) out of range, cannot parse \r\n $point"))
    elseif toreturn["Density(kg/m3)"] < 1000
        throw(ParsePerplexError("Density = $(toreturn["Density(kg/m3)"]) out of range, cannot parse \r\n $point"))
    end
    
    # Get temp and pressure 
    targets = ["P(bar)", "T(K)"]
    properties_i = findfirst(a->occursin("stable phases at", a), lowercase.(titles))
    rows = split(out[properties_i], r"[\r\n]+")
    for r in rows
        n = strip(split(r, "=")[1])
        v = parse(Float64, split(r, "=")[2])
        if !(n in targets)
            throw(ParsePerplexError("Got unexpected property $n, cannot parse \r\n $point"))
        else 
            toreturn[n] = v
        end
    end
    
    # Get Ks, mu, vp, vs, Ks_T(bar/K),  Ks_P,  Mu_T(bar/K),  Mu_P
    # These exist in separate sections
    targets = ["Ks(bar)", "Mu(bar)", "Vp(km/s)", "Vs(km/s)", "Ks_T(bar/K)",  "Ks_P",  "Mu_T(bar/K)",  "Mu_P"]
    for title in ["seismic properties", "seismic derivatives"]
        properties_i = findfirst(a->occursin(title, a), lowercase.(titles))
        if isnothing(properties_i)
            throw(ParsePerplexError("Could not find $title in output titles $titles. Did you forget to turn on seismic_output=all in perplex_options.dat?"))
        end
        rows = split(out[properties_i], r"[\r\n]+")
        header = split(rows[1], r"\s{2,}")
        rows = rows[2:end]
        cols = [strip(split(r, r"\s{2,}")[1]) for r in rows]
        if !("System" in cols)
            throw(ParsePerplexError("System not found in seismic properties, cannot parse point \r\n $point"))
        end
        systemi = findfirst(isequal("System"), cols)
        for t in targets
            if !(t in header)
                continue # may be in other title 
                #throw(ParsePerplexError("$t not found in seismic properties, cannot parse point \r\n $point"))
            end
            ti = findfirst(isequal(t), header)
            toreturn[t] = parse(Float64, split(rows[systemi], r"\s{2,}")[ti])
        end
    end 

    # Check that we've got em all 
    for target in targets
        if !(target in keys(toreturn))
            throw(ParsePerplexError("Did not find property $target"))
        end 
    end 
    
    return toreturn
    
end

# """
# For a perplex point output, return dictionary from endmember phases, including solutions, to the vol% of that phase 
# Unused, i think 
# """
# function parse_perplex_point(point, solutions)
#     if !solutions
#         return parse_perplex_point(point)
#     end
    
#     titles = [m.match for m in eachmatch(r"[\r\n]+[A-Za-z0-9\s)(/]+:[\r\n]+", point)]
#     out = split(point, r"[\r\n]+[A-Za-z0-9\s)(/]+:[\r\n]+")[2:end] # first is space above first title
#     toreturn = Dict()
    
#     # Parse phase compositions
#     phase_comp_i = findfirst(a->occursin("phase compositions", a), lowercase.(titles))
#     rows = split(out[phase_comp_i], r"[\r\n]+")
#     header = split(rows[1], r"\s{2,}")
#     names = [strip(split(r, r"\s{2,}")[1]) for r in rows[2:end]]
#     volp = zeros(length(names))
#     target = findfirst(isequal("vol %"), header)
#     for (i,r) in enumerate(rows[2:end])
#         volp[i] = parse(Float64, split(r, r"\s{2,}")[target])
#     end
    
#     # Add phases to toreturn. Some solutions may occur twice. 
#     for (i,e) in enumerate(names)
#         current = get(toreturn, e, 0)
#         toreturn[e] = current+ volp[i]
#     end
    
#     if abs(sum(values(toreturn)) - 100) > 1
#         throw(ParsePerplexError("Endmembers do not sum to 100"))
#     end

#     return toreturn
# end 









