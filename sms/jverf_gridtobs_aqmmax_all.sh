#! /bin/sh

set -x

STARTDATE=2020030700
ENDDATE=2020030700

DATE4=$STARTDATE

export DATE4

while [ $DATE4 -le $ENDDATE ]; do
#DATE=`cut -c 1-8 date4`

#export DATE

echo $DATE4 > date4

bsub < jverf_gridtobs_pmmax_06.sms
bsub < jverf_gridtobs_pmmax1_06.sms
sleep 200
bsub < jverf_gridtobs_aqmmax_12.sms
bsub < jverf_gridtobs_pmmax_12.sms
bsub < jverf_gridtobs_aqmmax1_12.sms
bsub < jverf_gridtobs_pmmax1_12.sms
sleep 500

DATEM=`/nwprod/util/exec/ndate -24 $DATE4`
echo $DATEM > datem
DATE=`cut -c 1-8 datem`
 
###rm -f -r /meso/save/wx20ps/com/verf/prod/gridtobs.$DATE

DATE4=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +24 $DATE4`
export DATE4

done

exit
