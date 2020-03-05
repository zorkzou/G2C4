# G2C4
G2C4 is an interface program between [Gaussian](http://www.gaussian.com/) and [CFour](http://www.cfour.de/) (ver. 2.00beta or 2.1) via the `external` keyword of [Gaussian](http://www.gaussian.com/). [MRCC](http://www.mrcc.hu/) may also be used through [CFour](http://www.cfour.de/).

## Compilation

    > F90 -O3 g2c4.f90 -o g2c4.exe

where `F90` can be `gfortran`, `pgf90`, `ifort`, or other Fortran 90 compilers.

## How to run CFour in Gaussian

1. Put g2c4.exe as well as 3 templet files scripts/cfour.templet-* in to a folder. For example, /home/myID/Chemsoft/G2C4
2. Put scripts/run-cfour.sh into your Gaussian calculation folder, and modify run-cfour.sh.
* Line 4 has to be modified if MRCC is going to be used,
* Lines 7 to 13 are related to CFour,
* Lines 15 to 28 are related to G2C4, where cfour_templet specifies the templet file used for CFour calculations. Change it to $g2c4dir/cfour.templet-ccsd, for example.
3. If necessary, modify the templet file, such as basis set, ECP, occupation, convergence criteria, and so on.
4. Make a copy of test/opt-freq.inp in your Gaussian calculation folder, and submit a Gaussian job.

