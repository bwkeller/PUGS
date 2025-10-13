#!/usr/bin/env bash

# Build our local venv
python3 -m venv .venv
source .venv/bin/activate
pip install -r ../pip-requirements.txt
pip install -r requirements.txt

# Run CAMB
camb inputs/planck_2018_CAMB.ini &> outputs/camb.log

# Run genetIC
cp inputs/genetIC_volume.txt .
docker run --rm -v `pwd`:/w/ --user $(id -u):$(id -g) apontzen/genetic:1.5.0 /w/genetIC_volume.txt


# Clean up
mv planck_2018_transfer_out.dat outputs/
mv tipsy.param outputs/
mv IC_output.params outputs/
rm planck_2018*
rm genetIC_volume.txt

