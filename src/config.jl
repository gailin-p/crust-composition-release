

# All composition elements 
COMPOSITION_ELEMENTS = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O_Total","CO2"]
export COMPOSITION_ELEMENTS
# Same elements, but perplex-ready format 
PERPLEX_COMPOSITION_ELTS = ["SIO2", "TIO2", "AL2O3", "FEO", "MGO", "CAO", "NA2O", "K2O", "H2O", "CO2"]
export PERPLEX_COMPOSITION_ELTS

# All elements resampled by resampleEarthChem
RESAMPLED_ELEMENTS = vcat(COMPOSITION_ELEMENTS, ["Latitude", "Longitude", "Age"])
export RESAMPLED_ELEMENTS

# All elements output by resampleEarthChem and sent to perplex 
PERPLEX_ELEMENTS = vcat(["index"], RESAMPLED_ELEMENTS, "geotherm", "upper", "middle", "lower")
export PERPLEX_ELEMENTS

UPPER="upper"
MIDDLE="middle"
LOWER="lower"
LAYER_NAMES = (UPPER, MIDDLE, LOWER)
export LAYER_NAMES

# ign file 
IGN_FILE = "igncn1.mat"