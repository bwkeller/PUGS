This package contains all the required scripts and input files used to generate the $50 h^{-1}$ Mpc, 
2048^3 collisionless volume that PUGS uses to generate zoom ICs from.

Directory Structure
-------------------
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
      conditions for the large, collisionless volume.  When run with genetIC,
      it should output a tipsy IC ready to run.
