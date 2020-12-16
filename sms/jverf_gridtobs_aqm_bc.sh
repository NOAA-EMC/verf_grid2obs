#! /bin/sh

set -x

STARTDATE=2020033100
ENDDATE=2020033100

DATE4=$STARTDATE

export DATE4

while [ $DATE4 -le $ENDDATE ]; do
#DATE=`cut -c 1-8 date4`

#export DATE

echo $DATE4 > date4

bsub < jverf_gridtobs_aqm1_00.sms
bsub < jverf_gridtobs_pm1_00.sms
sleep 100
bsub < jverf_gridtobs_aqm1_03.sms
bsub < jverf_gridtobs_pm1_03.sms
sleep 100
bsub < jverf_gridtobs_aqm1_06.sms
bsub < jverf_gridtobs_pm1_06.sms
sleep 100
bsub < jverf_gridtobs_aqm1_09.sms
bsub < jverf_gridtobs_pm1_09.sms
sleep 100
bsub < jverf_gridtobs_aqm1_12.sms
bsub < jverf_gridtobs_pm1_12.sms
sleep 100
bsub < jverf_gridtobs_aqm1_15.sms
bsub < jverf_gridtobs_pm1_15.sms
sleep 100
bsub < jverf_gridtobs_aqm1_18.sms
bsub < jverf_gridtobs_pm1_18.sms
sleep 100
bsub < jverf_gridtobs_aqm1_21.sms
bsub < jverf_gridtobs_pm1_21.sms
sleep 100
bsub < jverf_gridtobs_pmmax1_06.sms
sleep 200
bsub < jverf_gridtobs_aqmmax1_12.sms
bsub < jverf_gridtobs_pmmax1_12.sms
#sleep 500

DATEM=`/nwprod/util/exec/ndate -24 $DATE4`
echo $DATEM > datem
DATE=`cut -c 1-8 datem`
 
###rm -f -r /meso/save/wx20ps/com/verf/prod/gridtobs.$DATE

DATE4=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +24 $DATE4`
export DATE4

done

exit
