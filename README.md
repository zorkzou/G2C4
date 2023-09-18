# G2C4
G2C4 is an interface program between [Gaussian](http://www.gaussian.com/) and [CFour](http://www.cfour.de/) (ver. 2.00beta or 2.1) via the `External` keyword of [Gaussian](http://www.gaussian.com/). [MRCC](http://www.mrcc.hu/) may also be used through [CFour](http://www.cfour.de/).

## Recent Changes

09/19/2023

1. The script `run-cfour.sh` has been updated. The initial guess of Hartree-Fock may be read from an old `JOBARC` file by default to speed up SCF convergence.

04/21/2023

1. Bug fix for the `ONIOM` calculation by Gaussian 16.c.

11/18/2021

1. Reorder atoms back to the initial ordering if necessary.

2. Bug fix for the `ONIOM` EIn file by Gaussian 16.c.

02/25/2020

1. The first public version.

## Compilation

    > F90 -O3 g2c4.f90 -o g2c4.exe

where `F90` can be `gfortran`, `nvf90` (`pgf90`), `ifort`, or other Fortran 90 compilers.

## How to run CFour in Gaussian

1. Put g2c4.exe as well as 3 templet files scripts/cfour.templet-* into a folder. For example, /home/myID/Chemsoft/G2C4

2. Put scripts/run-cfour.sh into your Gaussian calculation folder, and modify run-cfour.sh.
* Line 4 has to be modified if MRCC is going to be used,
* Lines 7 to 13 are related to CFour,
* Lines 15 to 28 are related to G2C4, where cfour_templet specifies the templet file used for CFour calculations. Change it to $g2c4dir/cfour.templet-ccsd, for example.

3. Don't forget to add the execute permission to run-cfour.sh: chmod +x run-cfour.sh

4. If necessary, modify the templet file, such as basis set, ECP, occupation, convergence criteria, and so on. The basis set and ECP files may also be defined in the lines 30 to 32 in run-cfour.sh.

5. Make a copy of test/opt-freq.inp in your Gaussian calculation folder, modify the path of run-cfour.sh in `external=` if necessary, and submit it as a Gaussian job.

## Some tips

1. In the case of single layer system (i.e. without using the `ONIOM` keyword), the `External` keyword activates the geometry optimization procedure with MM microiterations, which converges quite slow. `OPT(zmat)` (for Z-matrix coordinates only!) and `OPT(nomicro)` converge much faster. The `ONIOM` calculation has no such a problem since redundand internal coordinates being used by default.

2. If both high and low layers in `ONIOM` are computed by CFour (or other external programs), Gaussian optimization procedure may report an error "MOMM, but no IMMCRS" in L103. You may use `OPT(nomicro)` or reset `IOP(1/18)` to 20 to solve the problem. For example, `# ONIOM(External='./run-cfour-high.sh':External='./run-cfour-low.sh') IOP(1/18=20) OPT`

## Known problems

1. In the case of `ONIOM` by Gaussian 16.a (and maybe earlier revisions), `OPT` and `FREQ` have to be calculated separately through this interface.

2. If CFour cannot read ECPDATA defined by `%ECPDATA=...`, this file has to be provided in run-cfour.sh (line 32).

