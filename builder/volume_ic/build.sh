#!/usr/bin/env bash

# Check arguments
TEST=0
NUM_THREADS=96
if [ "$#" -ne 0 ]; then
    if [ "$1" == "-t" ]; then
        TEST=1
        NUM_THREADS=4
    fi
    if [[ "$#" -ne 1 || "$1" == "-h" ]]; then
        echo "Usage: ./build.sh [OPTION]"
        echo "Build the initial conditions for the collisionless PUGS volume"
        echo ""
        echo "Arguments:"
        echo "  -h  print this help message"
        echo "  -t  build a smaller 128^3 volume for testing"
        exit
    fi
fi

mkdir -p outputs

# Run CAMB
camb ../../inputs/planck_2018_CAMB.ini &> outputs/camb.log

# Run genetIC
cp ../../inputs/genetIC_volume.txt .
if [ "$TEST" -eq 1 ]; then
    sed -i -e 's/2048/128/g' genetIC_volume.txt
fi
docker run -e OMP_NUM_THREADS=$NUM_THREADS --rm -v `pwd`:/w/ --user $(id -u):$(id -g) apontzen/genetic:1.5.0 /w/genetIC_volume.txt

# Fix the 32-bit tipsy header issue
if [ "$TEST" -eq 0 ]; then
    ./fix_header.py outputs/DM_volume.tipsy
fi

# Clean up
mv planck_2018_transfer_out.dat outputs/
mv tipsy.param outputs/
mv IC_output.params outputs/
rm planck_2018*
rm genetIC_volume.txt
