#!/usr/bin/env bash
# Build the PUGS Parquet halo catalog from a simulation's snapshots and AHF
# output.  Everything is derived directly from the AHF files and the snapshots;
# there is no intermediate database.
#
#   ./build_catalog.sh test_vars
set -e

# shellcheck source=/dev/null
source "$1"

: "${PUGS_SIMULATION:?set PUGS_SIMULATION to the directory holding the snapshots}"
: "${PUGS_CATALOG:?set PUGS_CATALOG to the output catalog directory}"

ARGS=()
if [ -n "$PUGS_MIN_PARTICLES" ]; then
	ARGS+=(--min-particles "$PUGS_MIN_PARTICLES")
fi
if [ -n "$PUGS_NAME" ]; then
	ARGS+=(--name "$PUGS_NAME")
fi

python -m pugs.export "$PUGS_SIMULATION" --output "$PUGS_CATALOG" "${ARGS[@]}"
