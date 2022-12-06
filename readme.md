## Usage

### Steps for basic run:

See test_end_to_end.sh for an example of a small run. 

src/resampleEarthChem.jl
src/runPerplex.jl
src/inversion_binned_geotherm.jl

## Organization

scripts/ holds the scripts for basic data processing.
src/ holds code for inversion models, seismic models, and data input
data/ (gitignored) stores data produced by script
resources/ holds data sources referenced by scripts/ and src/
visualization/ holds scripts to visualize data
