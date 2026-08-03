#!/usr/bin/env bash
set -e

# shellcheck source=/dev/null
source "$1"

RUNNER=(tangos)
MPI_FLAGS=()
if [ -n "$PUGS_NPROCS" ]
then
    RUNNER=(mpirun -np "$PUGS_NPROCS" tangos)
    MPI_FLAGS=(--backend mpi4py --load-mode server-shared-mem)
fi

tangos add --quicker --min-particles 500 "$SIM"

tangos import-properties hostHalo numSubStruct Mhalo Rhalo Rmax r2 mbp_offset com_offset Vmax v_esc sigV lambda lambdaE Lx Ly Lz b c Eax Eay Eaz Ebx Eby Ebz Ecx Ecy Ecz ovdens Ekin Epot SurfP Phi0 cNFW

tangos import-ahf-trees

tangos remove-duplicates

"${RUNNER[@]}" write finder_mass shrink_center max_radius R200 R500 M500 M200 --for "$SIM" --with-prerequisites "${MPI_FLAGS[@]}"

# write needed properties for only the last (z=0) snapshot
"${RUNNER[@]}" write zlib_ids Rvir_indices N_mm z_lmm z25_mass z50_mass z75_mass --with-prerequisites --for "$SIM" --latest "${MPI_FLAGS[@]}"
