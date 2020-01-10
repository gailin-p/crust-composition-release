# Test if the exclusion of F (fluid) from solution models fixes bad values

using StatGeochem
using Plots; gr();
using ArgParse
using DelimitedFiles

s = ArgParseSettings()
@add_arg_table s begin
    "--indices"
        help = "Indices of samples to run"
        arg_type = Int
        nargs = '*'
        required = true
end
parsed_args = parse_args(ARGS, s)

## Some mafic and felsic samples for testing fluid model
# 1503,53.83,1.08,16.08,9.8979,2.89,8.07,3.32,0.5,1.8386,0.4407,48 # mafic, high CAO
# 2,74.14,0.25,14.9,2.1146,0.5,4.8,3.15,0.08,1.8386,0.4407,32 # felsic
# 9,66.48,0.48,15.38,3.6442,2.33,3.11,3.49,3.89,1.8386,0.4407,41 # This one didn't have the initial issue
# 25,54.6,2.21,16.4,9.448,2.89,4.35,2.36,4.16,1.8386,0.4407,41 # mafic, low CAO
# 26,56.2,2.15,15.7,9.0881,2.68,6.06,2.44,2.78,1.8386,0.4407,41 # mafic, medium CAO
# 144,53.88,0.92,14.4,8.1253,6.93,8.49,2.88,1.11,1.8386,0.4407,42 # mafic, high CAO

# General options
perplexdir = "/Users/gailin/dartmouth/crustal_structure/Perple_X/"
scratchdir = "/Users/gailin/dartmouth/crustal_structure/perplexed_pasta/"
exclude = ""
dataset = "hpha02ver.dat"
elts = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]
P_range = [1, 20000]
dpdz = 2900. * 9.8 / 1E5 * 1E3
depth = 40.2 # mean 550C depth of all samples
geotherm = 550.0/depth/dpdz

# composition dataset
dat = readdlm("../ignmajors.csv",',')

# Solution options
solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\nF\n"
solutions_nof = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\n"
solutions_f_melt = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\nF\nmelt(HP)\n"


for index in parsed_args["indices"]
    sample = dat[index,:]
    comp = sample[2:11]


    perplex_configure_geotherm(perplexdir, scratchdir, comp, elements=elts,
        P_range=P_range, geotherm=geotherm, dataset=dataset,
        solution_phases=solutions_f_melt)

    modes = perplex_query_modes(perplexdir, scratchdir)

    # Where do we first see F?
    if "F" in modes["elements"]
        first = findfirst(!isnan, modes["F"])
        println("depth of first not-nan F mode $(modes["P(bar)"][first]/dpdz)")
    end

    # Modes is pretty long, make it smaller for graphing
    println(length(modes["T(K)"]))
    all = 1:1:length(modes["T(K)"])
    for elem in modes["elements"]
        modes[elem] = modes[elem][(all .% 10) .== 0]
    end
    println(length(modes["T(K)"]))

    # Plot averaged seismic properties
    depth = modes["P(bar)"] ./ dpdz

    linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    h = plot(xlabel="depth (km)", ylabel="Weight percent")
    for (l,m) in enumerate(modes["elements"][3:end])
        if sum(map(!isnan,modes[m])) > 0
            plot!(h, depth, modes[m], label=m, linestyle=linestyles[l%length(linestyles)+1])
        end
    end

    plot!(h, size=(800,800), ylims=(0,50), fg_color_legend=:white, framestyle=:box)
    savefig(h, "../output/low_seismic/$(string(index))_phase_modes_f-melt.pdf")
end
