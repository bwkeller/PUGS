#!/usr/bin/env bash

tangos add --min-particles 500 $1
tangos link --sims $1
tangos import-properties Mhalo Rhalo lambda r2 cNFW --for $1
tangos write finder_mass shrink_center --for $1
tangos write --latest --with-prerequisites z_half --for $1
