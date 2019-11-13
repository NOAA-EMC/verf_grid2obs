set -x

DATE=$1
domain=$2

echo $DATE > datem00

CYC=00

DATEM01=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +1 $DATE`
echo $DATEM01 > datem01
DATEM02=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +2 $DATE`
echo $DATEM02 > datem02
DATEM03=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +3 $DATE`
echo $DATEM03 > datem03
DATEM04=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +4 $DATE`
echo $DATEM04 > datem04
DATEM05=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +5 $DATE`
echo $DATEM05 > datem05
DATEM06=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +6 $DATE`
echo $DATEM06 > datem06

PDY00=`cut -c 1-8 datem00`
HH00=`cut -c 9-10 datem00`
PDY01=`cut -c 1-8 datem01`
HH01=`cut -c 9-10 datem01`
PDY02=`cut -c 1-8 datem02`
HH02=`cut -c 9-10 datem02`
PDY03=`cut -c 1-8 datem03`
HH03=`cut -c 9-10 datem03`
PDY04=`cut -c 1-8 datem04`
HH04=`cut -c 9-10 datem04`
PDY05=`cut -c 1-8 datem04`
HH05=`cut -c 9-10 datem05`
PDY06=`cut -c 1-8 datem06`
HH06=`cut -c 9-10 datem06`

DATE3=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate -3 $DATE`
echo $DATE3 > date3
PDY3=`cut -c 1-8 date3`
HH3=`cut -c 9-10 date3`


case $network in
  global)    bufrnet=meso;;
  globalak)  bufrnet=meso;;
  global40)  bufrnet=meso;;
  gdas)      bufrnet=meso;;
  gfse)      bufrnet=meso;;
  ozone)     bufrnet=aqm;;
  ozone1)    bufrnet=aqm;;
  aqm)       bufrnet=aqm;;
  aqm1)      bufrnet=aqm;;
  pm)        bufrnet=aqm;;
  pm1)       bufrnet=aqm;;
  meso)      bufrnet=meso;;
  satimg)    bufrnet=satimg;;
  hwrf)      bufrnet=global;;
   *)        bufrnet=meso;;
esac

if [ "$bufrnet" = "satimg" ]
then
  case $HH00 in
      00) ndas1=$COMG/gfs.${PDY00}/gfs.t${HH00}z.geoimr.tm00.bufr_d;;
      06) ndas1=$COMG/gfs.${PDY00}/gfs.t${HH00}z.geoimr.tm00.bufr_d;;
      12) ndas1=$COMG/gfs.${PDY00}/gfs.t${HH00}z.geoimr.tm00.bufr_d;;
      18) ndas1=$COMG/gfs.${PDY00}/gfs.t${HH00}z.geoimr.tm00.bufr_d;;
  esac

elif [ "$bufrnet" = "global" ]
then
  case $HH00 in
    00) ndas1=$COMG/gdas.$PDY00/gdas1.t${HH00}z.prepbufr
        ndas2=$COMG/$cyc/gfs.$PDY00/gfs.t${HH00}z.prepbufr
        ndas3=$COMN/ndas.$PDY12/ndas.t${HH12}z.prepbufr.tm12
        acars1=$COMN/ndas.$PDY12/ndas.t${HH12}z.prepbufr.acft_profiles_sfc.tm12;;
    03) ndas1=$COMN/ndas.$PDY09/ndas.t${HH09}z.prepbufr.tm09
        ndas2=$COMN/ndas.$PDY03/ndas.t${HH03}z.prepbufr.tm03;;
    06) ndas1=$COMG/gdas.$PDY00/gdas1.t${HH00}z.prepbufr
        ndas2=$COMG/gfs.$PDY00/gfs.t${HH00}z.prepbufr
        ndas3=$COMN/ndas.$PDY12/ndas.t${HH12}z.prepbufr.tm12;;
    09) ndas1=$COMN/ndas.$PDY09/ndas.t${HH09}z.prepbufr.tm09
        ndas2=$COMN/ndas.$PDY03/ndas.t${HH03}z.prepbufr.tm03;;
    12) ndas1=$COMG/gdas.$PDY00/gdas1.t${HH00}z.prepbufr
        ndas2=$COMG/gfs.$PDY00/gfs.t${HH00}z.prepbufr
        ndas3=$COMN/ndas.$PDY12/ndas.t${HH12}z.prepbufr.tm12;;
    15) ndas1=$COMN/ndas.$PDY09/ndas.t${HH09}z.prepbufr.tm09
        ndas2=$COMN/ndas.$PDY03/ndas.t${HH03}z.prepbufr.tm03;;
    18) ndas1=$COMG/gdas.$PDY00/gdas1.t${HH00}z.prepbufr
        ndas2=$COMG/gfs.$PDY00/gfs.t${HH00}z.prepbufr
        ndas3=$COMN/ndas.$PDY12/ndas.t${HH12}z.prepbufr.tm12;;
    21) ndas1=$COMN/ndas.$PDY09/ndas.t${HH09}z.prepbufr.tm09
        ndas2=$COMN/ndas.$PDY03/ndas.t${HH03}z.prepbufr.tm03;;
  esac
elif [ $bufrnet = "aqm" ]
then
  vday=`echo $DATE |cut -c1-8`
  vdaym1=`sh /gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.3/ush/finddate.sh $vday d-1`
  vdayp1=`sh /gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.3/ush/finddate.sh $vday d+1`

  case $domain in
    ozone) filnam=aqm.t12z.prepbufr.tm00
           COMBUFR=/gpfs/dell2/emc/verification/noscrub/Perry.Shafran/com/hourly/prod
           COMBUFR_IN1=$COMBUFR/hourly.$vdaym1
           COMBUFR_IN2=$COMBUFR/hourly.$vday;;
    ozone1) filnam=aqm.t12z.prepbufr.tm00
           COMBUFR=/gpfs/dell2/emc/verification/noscrub/Perry.Shafran/com/hourly/prod
           COMBUFR_IN1=$COMBUFR/hourly.$vdaym1
           COMBUFR_IN2=$COMBUFR/hourly.$vday;;
    pm)    filnam=aqm.t12z.anowpm.pb.tm024
           COMBUFR=/gpfs/dell2/emc/verification/noscrub/Perry.Shafran/com/hourly/prod
           COMBUFR_IN1=$COMBUFR/hourly.$vday
           COMBUFR_IN2=$COMBUFR/hourly.$vdayp1;;
    pm1)    filnam=aqm.t12z.anowpm.pb.tm024
           COMBUFR=/gpfs/dell2/emc/verification/noscrub/Perry.Shafran/com/hourly/prod
           COMBUFR_IN1=$COMBUFR/hourly.$vday
           COMBUFR_IN2=$COMBUFR/hourly.$vdayp1;;
  esac

  if [ $vcyc = "00" ]
  then
    ndas1=$COMBUFR_IN1/${filnam}
  else
    ndas1=$COMBUFR_IN2/${filnam} 
  fi

else
  case $HH00 in
    00) ndas1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.tm06
        ndas2=$COMN/nam.$PDY00/nam.t${HH00}z.prepbufr.tm00
        acars1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.acft_profiles_sfc.tm06;;
    01) ndas1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.tm05
        acars1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.acft_profiles_sfc.tm05;;
    02) ndas1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.tm04
        acars1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.acft_profiles_sfc.tm04;;
    03) ndas1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.tm03
        acars1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.acft_profiles_sfc.tm03;;
    04) ndas1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.tm02
        acars1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.acft_profiles_sfc.tm02;;
    05) ndas1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.tm01
        acars1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.acft_profiles_sfc.tm01;;
    06) ndas1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.tm06
        ndas2=$COMN/nam.$PDY00/nam.t${HH00}z.prepbufr.tm00
        acars1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.acft_profiles_sfc.tm06;;
    07) ndas1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.tm05
        acars1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.acft_profiles_sfc.tm05;;
    08) ndas1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.tm04
        acars1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.acft_profiles_sfc.tm04;;
    09) ndas1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.tm03
        acars1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.acft_profiles_sfc.tm03;;
    10) ndas1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.tm02
        acars1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.acft_profiles_sfc.tm02;;
    11) ndas1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.tm01
        acars1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.acft_profiles_sfc.tm01;;
    12) ndas1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.tm06
        ndas2=$COMN/nam.$PDY00/nam.t${HH00}z.prepbufr.tm00
        acars1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.acft_profiles_sfc.tm06;;
    13) ndas1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.tm05
        acars1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.acft_profiles_sfc.tm05;;
    14) ndas1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.tm04
        acars1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.acft_profiles_sfc.tm04;;
    15) ndas1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.tm03
        acars1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.acft_profiles_sfc.tm03;;
    16) ndas1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.tm02
        acars1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.acft_profiles_sfc.tm02;;
    17) ndas1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.tm01
        acars1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.acft_profiles_sfc.tm01;;
    18) ndas1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.tm06
        ndas2=$COMN/nam.$PDY00/nam.t${HH00}z.prepbufr.tm00
        acars1=$COMN/nam.$PDY06/nam.t${HH06}z.prepbufr.acft_profiles_sfc.tm06;;
    19) ndas1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.tm05
        acars1=$COMN/nam.$PDY05/nam.t${HH05}z.prepbufr.acft_profiles_sfc.tm05;;
    20) ndas1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.tm04
        acars1=$COMN/nam.$PDY04/nam.t${HH04}z.prepbufr.acft_profiles_sfc.tm04;;
    21) ndas1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.tm03
        acars1=$COMN/nam.$PDY03/nam.t${HH03}z.prepbufr.acft_profiles_sfc.tm03;;
    22) ndas1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.tm02
        acars1=$COMN/nam.$PDY02/nam.t${HH02}z.prepbufr.acft_profiles_sfc.tm02;;
    23) ndas1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.tm01
        acars1=$COMN/nam.$PDY01/nam.t${HH01}z.prepbufr.acft_profiles_sfc.tm01;;
  esac
fi

if [ -s $ndas1 ]
then
  cp $ndas1 prepda.${DATE}
  chmod 700 prepda.${DATE}
elif [ -s $ndas2 ]
  then
    cp $ndas2 prepda.${DATE}
    chmod 700 prepda.${DATE}
  else
  if [ -s $ndas3 ]
  then
    cp $ndas3 prepda.${DATE}
    chmod 700 prepda.${DATE}
  fi
fi

   if [ $NET != "aqm" -o $NET != "pm" ]
   then
     cat prepda.${DATE} $acars1 > prepda.temp
     mv prepda.temp prepda.${DATE}
   fi

ls -1 prepda.${DATE}
err3=$?

if [ $err3 -ne 0 ] ; then
  cp $COMG/gdas.${PDY00}/gdas1.t${HH00}z.prepbufr prepda.${DATE}
  err4=$?
  if [ $err4 -ne 0 ] ; then
    cp $COMG/gfs.${PDY00}/gfs.t${HH00}z.prepbufr prepda.${DATE}
  fi
fi

if [ $NET = "ngac" ]
then
cyc=$HH00                                   # center hour of dump (eg, '12' for 12z)
CYMD=${vday}$cyc
window_radius=1.5                            # 3 hour window == 1.5 hour radius
rc_dump=0
sh /gpfs/dell1/nco/ops/nwprod/obsproc_dump.v5.0.0/ush/dumpjb $CYMD $window_radius aodmod
rc_dump=$?
if [ $rc_dump -eq 0 ] ; then
  mv ../aodmod.ibm verf.t${cyc}z.aodmod.tm00.bufr_d
  mv ../aodmod.out verf.t${cyc}z.aodmod.tm00.bufr_d.out
fi                                           # rc_dump = 0

if [ -e verf.t${cyc}z.aodmod.tm00.bufr_d ]
then
 cp verf.t${cyc}z.aodmod.tm00.bufr_d prepda.${DATE}
else
 echo "FATAL ERROR: verf.t${cyc}z.aodmod.tm00.bufr_d NOT FOUND - Possible MODIS AOD data outage. Notify SDM and set this job to complete"
err=1
fi

fi

