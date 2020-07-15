#!/bin/ksh
#############################################################################
# Name of script:    verf_gridtobs_srefens.sh
# Purpose of script: This script processes the sref ensemble verification
# Arguments:
#  1. input parm card
#
# Notes: The following variables should be exported from the parent script
#
# vdate -- Verification Date, format yyyymmddhh
# vcyc  -- Verification Cycle, format hh
# frange -- Forecast Length
#
#############################################################################
set -x

if [ $# -lt 1 ]
then
   echo "Invalid argument"
   echo "usage: verf_gridtobs_srefens.sh $parmcard"
   err_exit
fi

input_parm=$1
domain=`echo $input_parm |awk -F"_" '{print $2}'`
group=`echo $input_parm |awk -F"_" '{print $3}'`
list=`echo $input_parm |awk -F"_" '{print $4}'`
seg=`echo $input_parm |awk -F"_" '{print $5}'`

# Execute the program:
export pgm=verf_gridtobs_gridtobs_ens
. prep_step

startmsg

mkdir -p $DATA/${domain}_${group}_${list}_${seg}
cp -d fort* $DATA/${domain}_${group}_${list}_${seg}
cp gridtobs_${domain}_${group}_${list}_${seg} $DATA/${domain}_${group}_${list}_${seg}
cd $DATA/${domain}_${group}_${list}_${seg}

ln -sf $PARMverf_gridtobs/verf_gridtobs.grid104 fort.20
ln -sf $PARMverf_gridtobs/verf_gridtobs.regions fort.21
ln -sf $DATA/${domain}_${group}_${list}_${seg}_${vdate}.vdb fort.50

$EXECverf_gridtobs/verf_gridtobs_gridtobs_ens <$DATA/$input_parm >gto.${domain}${vcyc}_${group}_${list}_${seg}.out
echo ${domain}_${group}_${list}_${seg}_${vdate};export err=$?; err_chk
grep "TOTAL" gto.${domain}${vcyc}_${group}_${list}_${seg}.out >>$COMOUT/srefens${XC}.t${vcyc}z.runlog

