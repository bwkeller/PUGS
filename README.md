# PUGS: The Portable Universal Galaxy Sampler [![build_volume_ic](https://github.com/bwkeller/PUGS/actions/workflows/build_volume_ic.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_volume_ic.yml) [![Catalog Build Test](https://github.com/bwkeller/PUGS/actions/workflows/build_catalog.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_catalog.yml) [![Container Build Test](https://github.com/bwkeller/PUGS/actions/workflows/build_container.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_container.yml) [![linter](https://github.com/bwkeller/PUGS/actions/workflows/linter.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/linter.yml) [![Deploy Docs](https://github.com/bwkeller/PUGS/actions/workflows/docs.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/docs.yml)

# Building the PUGS Data
## Requirements
- Python 3 (tested with version 3.9.16 and above)
- Docker (tested with version 28.5.1, build e180ab8)
## `inputs/
This is the directory containing the input param files for CAMB and genetIC
- `planck_2018_CAMB.ini`: this is the parameter file used to generate the
  CAMB transfer function used by genetIC to build the z=200 IC.  The default
  transfer function provided by genetIC is produced by a slightly different
  parameter file than the CAMB default, with the scalar spectra disabled
  (`get_scalar_cls = F`).  This produced a transfer function with a slightly
  too-low value of $\sigma_8$ (0.8117 vs. 0.8120 for the Plik best fit from
  Planck 2018).
- `genetIC_volume.txt`: parameter file used by genetIC to build the initial
  conditions for the large, collisionless volume.  When run with genetIC, it
  should output a tipsy IC ready to run.

## `builder/volume_ic`
This package contains all the required scripts and input files used to generate
the $50 h^{-1}$ Mpc, 2048^3 collisionless volume that PUGS uses to generate zoom
ICs from.  Be aware that because this IC is quite large, genetIC will require at
least 320 GB of memory to run.  It takes roughly 90 minutes to complete on 96
EPYC 9454 cores.

In order to get bitwise reproducibility, the IC _must_ be run with the same
number of threads.  If you are running this on a machine with fewer than 96
cores, the build will still succeed, but the performance may be poor.

### Directory Structure
- `build.sh`: This small bash script will run both CAMB and genetIC to build the
  PUGS volume.  Passing the -t argument will run in a testing mode, where a much
  smaller 128^3 volume will be built, and this build will use 4 OMP threads.
- `fix_header.py`: Because the PUGS volume contains 2^32 particles, it overflows
  the 32-bit integer that default tipsy format ICs use.  pynbody's tipsy reader
  does not currently support the sussheader-version introduced in pkdgrav3, so
  this script is needed to fix the header of the IC.
- `.checksums.txt`: MD5 checksums for the test volume.  This is used by the
  CI system to ensure that the test volume is bitwise-identical

## `builder/catalog`
This module builds the Parquet halo catalog from the snapshots produced by
running the volume IC through your N-body code, together with the AHF
halo-finder output beside them.  Every column is derived directly from the AHF
files and the snapshots; there is no intermediate database.

```bash
cd builder/catalog
source test_vars
./build_catalog.sh test_vars
```

The result is one Parquet file per snapshot plus a `provenance.json` sidecar,
readable with DuckDB, polars or pandas without installing PUGS.
