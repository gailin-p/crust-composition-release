

# All composition elements
const COMPOSITION_ELEMENTS = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2"]

# Same elements, but perplex-ready format
const PERPLEX_COMPOSITION_ELTS = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]

# Minor elements. Elts from R&G table 8 minus N, Ge, Br, Ru which do not exist in ign dataset.
const MINOR_ELTS = ["MnO", "P2O5", "Li", "Be", "B", "F", "S", "Cl", "Sc", "V", "Cr", "Co", "Ni", "Cu",
	"Zn", "Ga", "As", "Se", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Pd"]

REE = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]

# All elements resampled by resampleEarthChem
const RESAMPLED_ELEMENTS = vcat(COMPOSITION_ELEMENTS, MINOR_ELTS, REE, ["Latitude", "Longitude", "Age"])

# All elements output by resampleEarthChem and sent to perplex
const PERPLEX_ELEMENTS = vcat(["index"], RESAMPLED_ELEMENTS, "geotherm", "upper", "middle", "lower","exhumed", "formation_temp")

const UPPER="upper"
const MIDDLE="middle"
const LOWER="lower"
const LAYER_NAMES = (UPPER, MIDDLE, LOWER)

# ign file
const IGN_FILE = isfile("resources/igncn1.mat") ? "resources/igncn1.mat" : "../resources/igncn1.mat"

const dpdz = 2900. * 9.8 / 1E5 * 1E3; # bar/km.  density * gravity / 100,000  * 1000 m/km

##### Perple_X config
DEFAULT_DATASET = "hpha11ver.dat"
SOLUTIONS = "O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\nGlTrTsPg\nT\nB\nAnth\nChl(HP)"*
		"\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar"*
		"\nDo(HP)\n"
NPOINTS = 20
FLUID_ENDMEMBERS = "abL\nanL\ndiL\nenL\nfaL\nfliq\nfoL\nkspL\nmliq\nqL\nsiL\nq8L\nfa8L\nfo8L\nsil8L\nh2oL\nh2o8L\n"

DEFAULT_PERPLEX = "/Users/f0043n9/dartmouth/crustal_structure/perplex"
DEFAULT_SCRATCH = "/Users/f0043n9/dartmouth/crustal_structure/perplex_scratch"
