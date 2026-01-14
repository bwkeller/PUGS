#!/usr/bin/env bash
set -e

source config_vars

tangos add --min-particles 500 $SIM

tangos import-properties

tangos import-ahf-trees

tangos remove-duplicates

tangos write finder_mass shrink_center max_radius R500 M200 M500 --for $SIM --with-prerequisites

# write needed properties for only the last (z=0) snapshot
tangos write zlib_ids N_mm z_lmm z25_mass z50_mass z75_mass --with-prerequisites --for $SIM --latest
