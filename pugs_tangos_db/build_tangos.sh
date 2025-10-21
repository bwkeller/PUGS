#!/usr/bin/env bash

source .venv/bin/activate
source config_vars

tangos add --min-particles 500 $SIM
tangos link --sims $SIM
tangos import-properties Mhalo Rhalo lambda r2 cNFW --for $SIM
tangos write finder_mass shrink_center --for $SIM
