#!/usr/bin/env bash
# Build the PUGS Apptainer/Singularity container.
#
# The bundled halo catalog can be supplied with -c. If it is omitted, the
# catalog.placeholder directory is bundled instead so the container build still
# succeeds end-to-end before the 2048^3 catalog is ready.

set -euo pipefail

usage() {
	cat <<EOF
Usage: ./build_container.sh [-c CATALOG_DIR] [-o OUTPUT_SIF] [-h]

Options:
  -c CATALOG_DIR  Path to the PUGS Parquet halo catalog directory to bundle
                  into the container. Defaults to ./catalog.placeholder.
  -o OUTPUT_SIF   Output path for the built container.
                  Defaults to ./pugs.sif.
  -h              Show this help message and exit.
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CATALOG="${SCRIPT_DIR}/catalog.placeholder"
OUTPUT="${SCRIPT_DIR}/pugs.sif"

while getopts ":c:o:h" opt; do
	case "$opt" in
	c) CATALOG="$OPTARG" ;;
	o) OUTPUT="$OPTARG" ;;
	h)
		usage
		exit 0
		;;
	\?)
		echo "Unknown option: -$OPTARG" >&2
		usage >&2
		exit 2
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		usage >&2
		exit 2
		;;
	esac
done

if [ ! -d "$CATALOG" ]; then
	echo "ERROR: catalog directory not found: $CATALOG" >&2
	exit 1
fi

# Pick whichever container builder is on PATH. Apptainer is preferred on
# modern systems; Singularity is the historical name and still ships on
# many HPC sites.
if command -v apptainer >/dev/null 2>&1; then
	BUILDER=apptainer
elif command -v singularity >/dev/null 2>&1; then
	BUILDER=singularity
else
	echo "ERROR: neither apptainer nor singularity was found on PATH." >&2
	exit 1
fi

STAGED="${SCRIPT_DIR}/_catalog.staged"
cleanup() { rm -rf "$STAGED"; }
trap cleanup EXIT

rm -rf "$STAGED"
cp -R "$CATALOG" "$STAGED"

if [ "$CATALOG" = "${SCRIPT_DIR}/catalog.placeholder" ]; then
	echo "WARNING: building with the placeholder catalog. The bundled catalog" >&2
	echo "         holds no halo data. Re-run with -c /path/to/catalog once" >&2
	echo "         the 2048^3 catalog is ready." >&2
fi

echo "Builder: $BUILDER"
echo "Catalog: $CATALOG"
echo "Output:  $OUTPUT"

cd "$SCRIPT_DIR"
"$BUILDER" build --force "$OUTPUT" pugs.def
