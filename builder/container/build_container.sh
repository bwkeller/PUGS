#!/usr/bin/env bash
# Build the PUGS Apptainer/Singularity container.
#
# The bundled TANGOS database can be supplied with -d. If it is omitted, the
# pugs.db.placeholder file is bundled instead so the container build still
# succeeds end-to-end before the 2048^3 catalog is ready.

set -euo pipefail

usage() {
    cat <<EOF
Usage: ./build_container.sh [-d DB_PATH] [-o OUTPUT_SIF] [-h]

Options:
  -d DB_PATH      Path to the TANGOS SQLite database to bundle into the
                  container. Defaults to ./pugs.db.placeholder.
  -o OUTPUT_SIF   Output path for the built container.
                  Defaults to ./pugs.sif.
  -h              Show this help message and exit.
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_PATH="${SCRIPT_DIR}/pugs.db.placeholder"
OUTPUT="${SCRIPT_DIR}/pugs.sif"

while getopts ":d:o:h" opt; do
    case "$opt" in
        d) DB_PATH="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage >&2; exit 2 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage >&2; exit 2 ;;
    esac
done

if [ ! -f "$DB_PATH" ]; then
    echo "ERROR: database file not found: $DB_PATH" >&2
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

STAGED_DB="${SCRIPT_DIR}/_pugs.db.staged"
cleanup() { rm -f "$STAGED_DB"; }
trap cleanup EXIT

cp "$DB_PATH" "$STAGED_DB"

if [ "$DB_PATH" = "${SCRIPT_DIR}/pugs.db.placeholder" ]; then
    echo "WARNING: building with the placeholder DB. The bundled catalog will" >&2
    echo "         not be a valid SQLite file. Re-run with -d /path/to/pugs.db" >&2
    echo "         once the 2048^3 catalog is ready." >&2
fi

echo "Builder:  $BUILDER"
echo "Database: $DB_PATH"
echo "Output:   $OUTPUT"

cd "$SCRIPT_DIR"
"$BUILDER" build --force "$OUTPUT" pugs.def
