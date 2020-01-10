### Run perplex sample, plot seismic properties over depth
# UNFINISHED

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--ign"
        help = "ignmajors.csv perplex-ready data file"
        arg_type = String
        default = "../ignmajors.csv"
    "--indices"
        help = "Indices of samples to run"
        arg_type = Int
        nargs = '*'
        required = true
end
parsed_args = parse_args(ARGS, s)

# General options
perplexdir = "/Users/gailin/dartmouth/crustal_structure/Perple_X/"
scratchdir = "/Users/gailin/dartmouth/crustal_structure/perplexed_pasta/"
exclude = ""
dataset = "hpha02ver.dat"
elts = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]
P_range = [1, 20000]
dpdz = 2900. * 9.8 / 1E5 * 1E3

solutions = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\nF\n"
solutions_nof = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nDo(HP)\n"

# Plot all samples together
for i in parsed_args["indices"]
    # Run sample

    # For now, just run vp, could add command line options in the future
end
