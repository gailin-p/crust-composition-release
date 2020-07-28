

# All composition elements 
const COMPOSITION_ELEMENTS = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2"]

# Same elements, but perplex-ready format 
const PERPLEX_COMPOSITION_ELTS = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]

# All elements resampled by resampleEarthChem
const RESAMPLED_ELEMENTS = vcat(COMPOSITION_ELEMENTS, ["Latitude", "Longitude", "Age"])

# All elements output by resampleEarthChem and sent to perplex 
const PERPLEX_ELEMENTS = vcat(["index"], RESAMPLED_ELEMENTS, "geotherm", "upper", "middle", "lower","exhumed")

const UPPER="upper"
const MIDDLE="middle"
const LOWER="lower"
const LAYER_NAMES = (UPPER, MIDDLE, LOWER)

# ign file 
const IGN_FILE = "igncn1.mat"

