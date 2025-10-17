This package contains all the required scripts and input files used to generate
the $50 h^{-1}$ Mpc, 2048^3 collisionless volume that PUGS uses to generate zoom
ICs from.  Be aware that because this IC is quite large, genetIC will require at
least 320 GB of memory to run.  It takes roughly 90 minutes to complete on 96
EPYC 9454 cores.

Requirements
------------
- Python 3 (tested with version 3.9.16 and above)
- Docker (tested with version 28.5.1, build e180ab8)

Directory Structure
-------------------
- `requirements.txt`: python libraries (pinned to specific versions for
  reproducibility) required for this module of PUGS
-`inputs/`: This is the directory containing the input param files for CAMB and
  genetIC
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
- `build.sh`: This small bash script will run both CAMB and genetIC to build the
  PUGS volume.  Passing the -t argument will run in a testing mode, where a much
  smaller 128^3 volume will be built
- `fix_header.py`: Because the PUGS volume contains 2^32 particles, it overflows
  the 32-bit integer that default tipsy format ICs use.  TANGOS does not
  currently support the sussheader-version introduced in pkdgrav3, so this
  script is needed to fix the header of the IC.
