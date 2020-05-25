## Usage

Specify a data folder for the run:
1. Scripts will look for required data files here
2. All intermediate data products will go here.

### Steps for basic run:
src/resampleEarthChem.jl 
src/runPerplex.jl
src/inversion_binned_geotherm.jl 

## Organization

src/ holds the scripts for basic data processing.
Each script outputs an intermediate data file.
Each script requires a data folder with the output data file from the previous script.

visualization/ holds scripts to visualize data



