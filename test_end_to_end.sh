#!/bin/bash
set -e # Fail if any command fails

# Run an end-to-end test of base inversion model with very small N
# (100 composition samples, 1000 points to invert, 2 geotherm bins)

# REQUIRES:
#  - Perplex installed
#  - Change location of Perplex dir below to where your perplex is installed
#  - StatGeochem julia library

rm -rf data/test

# Resample composition files
julia scripts/resampleEarthChem.jl -o test -n 100 -w latlongage #--exhume true -b 2

# Perplex for each geotherm bin
mpiexec -np 3 --oversubscribe julia scripts/runPerplex.jl -d test -b 1 --scratch "/Users/f0043n9/dartmouth/crustal_structure/perplex_scratch" --perplex "/Users/f0043n9/dartmouth/crustal_structure/perplex" --perplex_dataset hpha11ver.dat 
#mpiexec -np 3 --oversubscribe julia scripts/runPerplex.jl -d test -b 2 --scratch "/Users/gailin/dartmouth/crustal_structure/perplexed_pasta/scratch" --perplex "/Users/gailin/resources/perplex-stable/"

# Invert
julia scripts/inversion_binned_geotherm.jl -d test -n 3000 --name default
julia visualization/age_comp.jl -p data/test/default
julia visualization/visualize_range_model.jl -d test --name default
