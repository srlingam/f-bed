#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - Feb. 2011
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
#if [ -d "$casePath/CFD/constant/polyMesh/boundary" ]; then
#    echo "mesh was built before - using old mesh"
#else
#    echo "mesh needs to be built"
#    cd $casePath/CFD
#    blockMesh
#fi

#- run parallel CFD-DEM in new terminal
gnome-terminal --title='cfdemSolverPiso settlingTest CFD'  -e "bash $casePath/parCFDDEMrun.sh" 
