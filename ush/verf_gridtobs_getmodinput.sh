#!/bin/ksh
#############################################################################
# Name of script:    verf_gridtobs_getmodinput.sh
# Purpose of script: This script obtains a list of model forecast output
#                     files to be fitted 
# Arguments:
#  1. domain: the domain for a  specific model, e.g. nam,nam40,namak,nam40ak
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
   echo "usage: verf_gridtobs_getmodinput.sh $region"
   err_exit
fi

wgrib=${wgrib:-$WGRIB}
copygb=${copygb:-$COPYGB}
wgrib2=${wgrib2:-$WGRIB2}
cnvgrib=${cnvgrib:-$CNVGRIB}
copygb2=${copygb2:-$COPYGB2}


domain=$1

cp $PARMverf_gridtobs/verf_gfs_parmlist gfs_parmlist

rm -rf input_$domain hour_$domain GRD* AWIP3D* nonzero_$domain

#####################################
# Determine the first available file
#####################################
echo "getmodinput" $vcyc ${domain}.cycles
grep "$vcyc" ${domain}.cycles
err_grep=$?

## If the vcyc is not one of the model run cycles, go back to the most recent cycle to 
## get the first input file
if [ $err_grep -ne 0 ]
then
   ifhr=192
   for pcyc in `cat ${domain}.cycles`
   do
     pfhr=`expr $vcyc - $pcyc`
     echo $pcyc $pfhr
     if [ $pfhr -lt 0 ]; then pfhr=`expr $pfhr + 24`; fi
     if [ $pfhr -le $ifhr ]; then ifhr=$pfhr; fi
   done
else
   ifhr=00
fi

## For some model, the forecast starting point is not zero, such as for the DGEX model, sfhr=84
fhr=`expr $ifhr + $sfhr`
fhra=$fhr
if [ $fhr -lt 10 ]
then
fhr=0$fhr
fhra=:$fhra
fi
export fhr

itr=1
ifire=0
while [ $fhr -le $frange ]
do
  if [ $model = ndas ]
  then
    adate=`$NDATE +$ifhr $vdate`
  else
    adate=`$NDATE -$fhr $vdate`
  fi
  aday=`echo $adate |cut -c1-8`
  acyc=`echo $adate |cut -c9-10`

  grep "$acyc" ${domain}.cycles
  err_grep=$?
  
  if [ $err_grep -eq 0 ]
  then
    echo $model $domain
    if [ $model = "nam" -o $model = "rap" -o $model = "conusnest" -o $model = "hrrr" -o $model = "hiresw" -o $model = "aknest" -o $model = "nssl4arw" ]
    then
    if [ $fhr -le 36 ]
    then
      if [ -e ${COMN}/nam.${aday}/nam.t${acyc}z.firewxnest.hiresf${fhr}.tm00.grib2 ]
      then
       cp ${COMN}/nam.${aday}/nam.t${acyc}z.firewxnest.hiresf${fhr}.tm00.grib2 firewx${fhr}.grb
       if [ -e firewx${fhr}.grb ]
       then
         let "ifire=ifire+1"
###         ${GRBINDEX} firewx${fhr}.grb firewx${fhr}i.grb
         cp ${COMN}/nam.${aday}/nam.t${acyc}z.firewxnest.hiresf${fhr}.tm00.grib2.idx firewx${fhr}i.grb
      fi
      fi
      fi
      fi
    case $model in
      fwis) cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr}
      ;;
      nam) cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr}
           pwd
           ls -l $DATA/${vcyc}
           ;;
      conusnest)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.conusnest.${filnam1}${fhr}${tmkk} GRD${fhr}
           ;;
      fv3nest) cp ${DIRIN}.${aday}/${acyc}/${runnam}.t${acyc}z.${filnam1}.f${fhr}.grib2 GRD${fhr}
           ;;
      aknest)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.alaskanest.${filnam1}${fhr}${tmkk} GRD${fhr}
           ;;
      hawaiinest)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.hawaiinest.${filnam1}${fhr}${tmkk} GRD${fhr}
                   gr="0 6 0 0 0 0 0 0 322 202 0 0 18000000 198000000 48 23025000 206025000 25000 25000 64"
                   $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
                   cp GRD${fhr}.1 GRD${fhr}
           ;;
      priconest)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.priconest.${filnam1}${fhr}${tmkk} GRD${fhr}
                   gr="0 6 0 0 0 0 0 0 171 121 0 0 16750000 291750000 48 19750000 296000000 25000 25000 64"
                   $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
                   cp GRD${fhr}.1 GRD${fhr}
           ;;
      ndas)fhour=`expr 12 - $ifhr`
           if [ $fhour -lt 10 ]; then fhour=0$fhour; fi
           aday=`$NDATE +$ifhr $vdate |cut -c1-8`
           acyc=`$NDATE +$ifhr $vdate |cut -c9-10`
           aday_p=`$NDATE +$fhour $vdate |cut -c1-8`
           acyc_p=`$NDATE +$fhour $vdate |cut -c9-10`

           if [ $itr -eq 1 ]
           then
             # Set fhr to 0 for the first iteration
             fhr=00
             cp $DIRIN.${aday_p}/${runnam}.t${acyc_p}z.${filnam1}${fhr}.tm${fhour}${tmkk} GRD${fhr}
           else
             if [ $ifhr -eq 00 ]
             then
               cp $DIRIN.${aday}/${runnam}.t${acyc}z.${filnam1}${fhr}.tm03${tmkk} GRD${fhr}
             else
               cp $DIRIN.${aday_p}/${runnam}.t${acyc_p}z.${filnam1}${fhr}.tm12${tmkk} GRD${fhr}
             fi
           fi
           ;;
      rap) cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1}f${fhr}${tmkk} GRD${fhr} 
           ;;
      narre) cp ${DIRIN}.${aday}/ensprod/${runnam}.t${acyc}z.${filnam1}.f${fhr}${tmkk} GRD${fhr}
           ;;
      rtma)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1} GRD${fhr}
           ;;
      rtma2p5)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1} GRD${fhr}
           ;;
      urma2p5)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1} GRD${fhr}
           ;;
      akrtma)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1} GRD${fhr}
           ;;
      hirtma)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1} GRD${fhr}
             gr="0 6 0 0 0 0 0 0 322 202 0 0 18000000 198000000 48 23025000 206025000 25000 25000 64"
             $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
             cp GRD${fhr}.1 GRD${fhr}
           ;;
      prrtma)cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1} GRD${fhr}
              gr="0 6 0 0 0 0 0 0 171 121 0 0 16750000 291750000 48 19750000 296000000 25000 25000 64"
              $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
             cp GRD${fhr}.1 GRD${fhr}
           ;;
      smartinit) cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr} 
                   if [ $domain = "smarthi" ]
                   then
                     gr="0 6 0 0 0 0 0 0 322 202 0 0 18000000 198000000 48 23025000 206025000 25000 25000 64"
                     $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
                     cp GRD${fhr}.1 GRD${fhr}
                 elif [ $domain = "smartpr" ]
                 then
                   gr="0 6 0 0 0 0 0 0 171 121 0 0 16750000 291750000 48 19750000 296000000 25000 25000 64"
                   $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
                   cp GRD${fhr}.1 GRD${fhr}
                 fi 
            ;;       
      sref) if [ $domain = "srmean" ]
            then
               if [ $fhr -eq 00 ]
               then
                 gribfile=${DIRIN}.${aday}/${acyc}/ensprod/${runnam}.t${acyc}z.${filnam1}.grib2
                 $WGRIB2 $gribfile | grep "anl" | $WGRIB2 -i $gribfile -grib GRD${fhr}.1
                 cp GRD${fhr}.1 GRD${fhr}
               else
                 gribfile=${DIRIN}.${aday}/${acyc}/ensprod/${runnam}.t${acyc}z.${filnam1}.grib2
		 echo "srmean" $gribfile
                 $WGRIB2 $gribfile | grep "${fhra} hour fcst" | $WGRIB2 -i $gribfile -grib GRD${fhr}.1
                 cp GRD${fhr}.1 GRD${fhr}
               fi
             else
              if [ -e ${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f${fhr}.grib2 ]
              then
               cp ${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f${fhr}.grib2 GRD${fhr}
              else
               cp ${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f${fhr} GRD${fhr}
              fi
             fi
            ;;
      srefak) if [ $domain = "srmeanak" ]
            then
               if [ $fhr -eq 00 ]
               then
                 gribfile=${DIRIN}.${aday}/${acyc}/ensprod/${runnam}.t${acyc}z.${filnam1}.grib2
                 $WGRIB2 $gribfile | grep "anl" | $WGRIB2 -i $gribfile -grib GRD${fhr}.1
                 cp GRD${fhr}.1 GRD${fhr}
               else
                 gribfile=${DIRIN}.${aday}/${acyc}/ensprod/${runnam}.t${acyc}z.${filnam1}.grib2
                 $WGRIB2 $gribfile | grep "${fhra} hour fcst" | $WGRIB2 -i $gribfile -grib GRD${fhr}.1
                 cp GRD${fhr}.1 GRD${fhr}
               fi
             else
              if [ -e ${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f${fhr}.grib2 ]
              then
               cp ${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f${fhr}.grib2 GRD${fhr}
              else
               cp ${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f${fhr} GRD${fhr}
              fi
             fi
            ;;
      srefx) if [ $domain = "srmeanx" ]
             then
               if [ $fhr -eq 00 ]
               then
                gribfile=${DIRIN}.${aday}/${acyc}/ensprod/${runnam}.t${acyc}z.${filnam1}.grib2
                $WGRIB2 $gribfile | grep "anl" | $WGRIB2 -i $gribfile -grib GRD${fhr}
               else
                gribfile=${DIRIN}.${aday}/${acyc}/ensprod_biasc/${runnam}.t${acyc}z.${filnam1}.grib2
                $WGRIB2 $gribfile | grep "${fhra} hour fcst" | $WGRIB2 -i $gribfile -grib GRD${fhr}
               fi
             else
                if [ $fhr -eq 00 ]
                then
                 gribfile=${DIRIN}.${aday}/${acyc}/pgrb/${runnam}.t${acyc}z.${filnam1}.f00.grib2
                 $WGRIB2 $gribfile | grep "anl" | $WGRIB2 -i $gribfile -grib GRD${fhr}
                else
                 gribfile=${DIRIN}.${aday}/${acyc}/pgrb_biasc/${runnam}.t${acyc}z.${filnam1}.grib2
                 $WGRIB2 $gribfile | grep "${fhra} hour fcst" | $WGRIB2 -i $gribfile -grib GRD${fhr}
                fi
            fi
           ;;
       ngac) fhrg=$(printf %03d ${fhr#0})
             cp ${DIRIN}.${aday}/${acyc}/chem/pgrb2ap25/gefs.chem.t${acyc}z.${filnam1}.f${fhrg}.grib2 GRD${fhr}
           ;;
      gfs) if [ $domain = "gfs" -a $GRIB = "grib2" ]
           then
             fhrg=$(printf %03d ${fhr#0})
           else
            # fhrg=$fhr
            fhrg=$(printf %03d ${fhr#0})
           fi
           if [ -e ${DIRIN}.${aday}/${acyc}/atmos/${runnam}.t${acyc}z.${filnam1}${fhrg}${tmkk} ]
           then
             cp ${DIRIN}.${aday}/${acyc}/atmos/${runnam}.t${acyc}z.${filnam1}${fhrg}${tmkk} TMP${fhr}
             wgrib2 TMP${fhr} | grep  -F -f gfs_parmlist | $WGRIB2 -i -grib  GRD${fhr} TMP${fhr}
             rm -f -r TMP${fhr}
             if [ $domain = "gfs12" -o $domain = "gfs12ak" ]
             then
             wgrib2 GRD${fhr} -set_grib_type same -new_grid_winds earth -new_grid ${cgrid} GRD${fhr}.1
             mv GRD${fhr}.1 GRD${fhr}
             fi
           fi
           ;;
     gfse) 
           fhrg=$(printf %03d ${fhr#0}) 
           if [ -e ${DIRIN}.${aday}/${acyc}/${runnam}.t${acyc}z.${filnam1}${fhrg}${tmkk} ]
           then
             # cp ${DIRIN}.${aday}/${acyc}/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} TMP${fhr}
             cp ${DIRIN}.${aday}/${acyc}/${runnam}.t${acyc}z.${filnam1}${fhrg}${tmkk} TMP${fhr}
             $WGRIB2 TMP${fhr} | grep  -F -f gfs_parmlist | $WGRIB2 -i -grib  GRD${fhr} TMP${fhr}
             rm -f -r TMP${fhr}
             $WGRIB2 GRD${fhr} -set_grib_type same -new_grid_winds earth -new_grid ${cgrid} GRD${fhr}.1
             mv GRD${fhr}.1 GRD${fhr}
            fi
             ;;
     gdas) if [ $domain = "gdas" -a $GRIB = "grib2" ]
           then
           fhrg=$(printf %03d ${fhr#0})
           else
            fhrg=$fhr
           fi
           if [ -e ${DIRIN}.${aday}/${acyc}/atmos/${runnam}.t${acyc}z.${filnam1}${fhrg}${tmkk} ]
           then
             cp ${DIRIN}.${aday}/${acyc}/atmos/${runnam}.t${acyc}z.${filnam1}${fhrg}${tmkk} GRD${fhr}
           fi
           ;;
      nssl4arw) if [ -e cp $DCOMROOT/prod/${aday}/wgrbbul/nssl_wrf/${filnam1}_${aday}${acyc}.f${fhr} ]
                then
                 cp $DCOMROOT/prod/${aday}/wgrbbul/nssl_wrf/${filnam1}_${aday}${acyc}.f${fhr} GRD${fhr}
                 if [ -e anl.grd ]
                 then
                  rm anl.grd
                 else
                  cp $DCOMROOT/prod/${aday}/wgrbbul/nssl_wrf/${filnam1}_${aday}${acyc}.f${fhr} anl.grd
                  if [ -e anl.grd ]
                  then
                   $CNVGRIB -g12 anl.grd anl.grd_grib2
                   mv anl.grd_grib2 anl.grd
                  fi
                 fi
                fi
           ;;
         hiresw) cp $DIRIN.${aday}/hiresw.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr}
             if [ $domain = "hiarw2" -o $domain = "hiarw" -o $domain = "hifv3" ]
               then
                gr="0 6 0 0 0 0 0 0 322 202 0 0 18000000 198000000 48 23025000 206025000 25000 25000 64"
                $COPYGB2 -g"$gr" -x GRD${fhr} GRD${fhr}.1
                mv GRD${fhr}.1 GRD${fhr}
             fi
            ;;
      hrrr) if [ $domain = hrrrak ]
             then
             cp ${DIRIN}.${aday}/alaska/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr}
             else
             cp ${DIRIN}.${aday}/conus/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr}
             fi
           ;;
         pm1) cp $DIRIN.${aday}/${runnam}.t${acyc}z.${filnam1}.f${fhr}.148.grib2 GRD${fhr}
           ;;
         aqm1) cp $DIRIN.${aday}/${runnam}.t${acyc}z.${filnam1}.f${fhr}.148.grib2 GRD${fhr}
           ;;
         pm) cp $DIRIN.${aday}/${runnam}.t${acyc}z.${filnam1}.f${fhr}.148.grib2 GRD${fhr}
           ;;
         aqm) cp $DIRIN.${aday}/${runnam}.t${acyc}z.${filnam1}.f${fhr}.148.grib2 GRD${fhr}
           ;;
        *)    if [ -e ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} ]
               echo "here"
               then 
               cp ${DIRIN}.${aday}/${runnam}.t${acyc}z.${filnam1}${fhr}${tmkk} GRD${fhr}
              fi
           ;;
    esac

    if [ $domain = "nssl4arw" ]
    then
     $CNVGRIB -g12 GRD${fhr} GRD${fhr}_grib2
     mv GRD${fhr}_grib2 GRD${fhr}
    fi

    if [ "$DOCYGB" = "YES" ]
    then
     if [ $GRIB = "grib1" ]
     then
      ${GRBINDEX} GRD${fhr} GRD${fhr}i
      $utilexec/copygb -g$cgrid GRD${fhr} GRD${fhr}i AWIP3D${fhr}.tm00
     fi
     if [ $GRIB = "grib2" ]
     then
      $GRB2INDEX GRD${fhr} GRD${fhr}i
      $COPYGB2 -g"$cgrid" GRD${fhr} GRD${fhr}i AWIP3D${fhr}.tm00
#      ${COPYGB2} -g"$cgrid" -x GRD${fhr} AWIP3D${fhr}.tm00
     fi
    else
      if [ -e GRD${fhr} ]
        then
        cp GRD${fhr} AWIP3D${fhr}.tm00
      fi
    fi

    if [ -e AWIP3D${fhr}.tm00 ] 
    then
     if [ $GRIB = "grib2" ]
     then
      $GRB2INDEX AWIP3D${fhr}.tm00 AWIP3D${fhr}i.tm00
     fi
     if [ $GRIB = "grib1" ]
     then
      $GRBINDEX AWIP3D${fhr}.tm00 AWIP3D${fhr}i.tm00
     fi
    fi

#    if [ $fhr -eq 00 -o $fhr -eq 03 ]

    if [ $model = "nssl4arw" ]
    then

     if [ -e anl.grd -a $itr -eq 1 ]
     then
      cp anl.grd ANL.tm00
      $GRB2INDEX ANL.tm00 ANLi.tm00

     fi

    else

    if [ $itr -eq 1 ]
    then

    cp AWIP3D${fhr}.tm00 ANL.tm00
    cp AWIP3D${fhr}i.tm00 ANLi.tm00

    fi

    fi
    
###  For fire weather  ###
 
    if [ -e firewx${fhr}.grb ]
    then

    if [ $fhr -le 36 ]
    then
    echo "firewx${fhr}.grb" > firewx
    echo "firewx${fhr}i.grb" >> firewx
    echo "maskout.grb" >> firewx

    export XLFUNIT_30=firewx${fhr}.grb
    export XLFUNIT_31=firewx${fhr}i.grb
    export XLFUNIT_32=maskout.grb

    startmsg
    $EXECverf_gridtobs/verf_gridtobs_writemask < firewx >>writemask${fhr}.out
    export err=$?; err_chk

    $GRB2INDEX maskout.grb maskouti.grb
    cp maskout.grb firemask${ifire}.grb
    cp maskouti.grb firemask${ifire}i.grb

    echo $fhr >>hour_${domain}_fwis

    fi
    fi

    fhri=$fhr
    if [ $fhr -lt 100 ]; then fhri=0$fhr; fi

    # generate order list:
    if [ -s AWIP3D${fhr}.tm00 ]
    then
     if [ $model = "ngac" ]
      then
       echo "AWIP3D${fhr}.tm00" >>nonzero_${domain}
       echo "$fhri $pll3 AWIP3D${fhr}.tm00" >>input_${domain}
       echo "$fhri $pll3 AWIP3D${fhr}i.tm00" >>input_${domain}
       echo $fhr >>hour_${domain}
      else
       echo "AWIP3D${fhr}.tm00" >>nonzero_${domain}
       echo "$fhri $pll3 AWIP3D${fhr}.tm00" >>input_${domain}
       echo "$fhri $pll3 AWIP3D${fhr}i.tm00" >>input_${domain}
       echo $fhr >>hour_${domain}
     fi
    fi
  fi

  if [ "$fhr_next" -gt 0 ]
  then 
     fhr=$fhr_next
     fhra=$fhr
  else
     let "fhr=fhr+$inc"
     fhra=$fhr
  fi
  let "itr=itr+1"
  if [ $fhr -lt 10 ]
  then
   fhr=0$fhr
   fhra=:$fhra
  fi

done
