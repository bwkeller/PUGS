#!/usr/bin/env bash

source .venv/bin/activate
source config_vars

tangos add --min-particles 500 $SIM

tangos write max_radius shrink_center --for $SIM


# write needed properties for the last (z=0) snapshot
tangos write zlib_ids --for $SIM --latest

#tangos link --sims $SIM
#tangos import-properties Mhalo Rhalo lambda r2 cNFW --for $SIM
#tangos write finder_mass shrink_center --for $SIM
