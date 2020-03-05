#!/bin/bash

# MRCC (optional)
#export PATH=/home/Chemsoft/MRCC/bin:$PATH

# CFour 2
export PATH=/home/Chemsoft/CFour_2.1/bin:$PATH
export C4SCR=/tmp/Scratch-CFour
export OMP_NUM_THREADS=2
#
rm -rf $C4SCR
mkdir $C4SCR
cd $C4SCR

# G2C4
export g2c4dir=/home/Chemsoft/G2C4
export cfour_scr=$C4SCR
#
export cfour_templet=$g2c4dir/cfour.templet
export gaussian_ein=$2
export gaussian_eou=$3
export cfour_zmt=$cfour_scr/ZMAT
export cfour_out=$cfour_scr/test.out
export cfour_grd=$cfour_scr/GRD
export cfour_dip=$cfour_scr/DIPOL
export cfour_fcm=$cfour_scr/FCMFINAL
export cfour_apt=$cfour_scr/DIPDER
export cfour_pol=$cfour_scr/POLAR

# generate CFour input file
$g2c4dir/g2c4.exe -g2c -gin $gaussian_ein -ctp $cfour_templet -cin $cfour_zmt

xcfour > $cfour_out

# write Gaussian's *.EOu file
$g2c4dir/g2c4.exe -c2g -gin $gaussian_ein -gou $gaussian_eou -cou $cfour_out \
  -c_grd $cfour_grd -c_dip $cfour_dip \
  -c_fcm $cfour_fcm -c_apt $cfour_apt -c_pol $cfour_pol


