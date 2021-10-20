#!/bin/ksh
######################################################################
#  UNIX Script Documentation Block
#                      .
# Script name:         verif_gridtobs_aqmmax.sh
# Script description:  Run max (04z-04z) ozone verification for CMAQ
#
# Author:      Marina Tsidulko     Org: NP22       Date: 2006-02-08
#
# Abstract: This script runs verification for max ozone concentrations
#           for CMAQ both for 1hr max and 8hr max
#
# Script history log:
# 2006-02-08    Marina Tsidulko
# 2008-08-04    Julia Zhu
#               Modification for production implementation
######################################################################
set -xa
msg="JOB $job HAS BEGUN"
postmsg "$jlogfile" "$msg"

model=$1
bc=$2

wgrib2=${wgrib2:-$WGRIB2}

if [ $bc -eq 0 ]
then
model_p="cmaq5x"
cmaq_domain="cmaq5x"
bctag=""
elif [ $bc -eq 1 ]
then
model_p="cmaq5xbc"
cmaq_domain="cmaq5xbc"
bctag="_bc"
fi
#### To process the Ozone part only--Julia 20090514
#regions="ozonemax pmmax pmave"
if [ $model = "aqmmax" -o $model = "aqmmax1" ] 
then
regions="ozonemax"
elif [ $model = "pmmax" -o $model = "pmmax1" ]
then
regions="pmave pmmax"
fi
#####regions="ozonemax pmave pmmax"

vdaym1=`$NDATE -24 $vdate |cut -c1-8`
vdaym2=`$NDATE -48 $vdate |cut -c1-8`
vdayp1=`$NDATE +24 $vdate |cut -c1-8`
vdayp2=`$NDATE +48 $vdate |cut -c1-8`
vdatep1=`$NDATE +24 $vdate`

if [ $SENDCOM = YES ]
then
  if [ ! -s $COMVSDB/${model_p} ]; then mkdir -p $COMVSDB/${model_p}; fi
#  if [ -s $COMVSDB/${model_p}/${cmaq_domain}max_${vday}.vsdb ]
#  then
#    sed -e "/$vday$vcyc/d" $COMVSDB/${model_p}/${cmaq_domain}max_${vday}.vsdb >$COMVSDB/${model_p}/${cmaq_domain}max_${vday}.vsdb1
#    mv $COMVSDB/${model_p}/${cmaq_domain}max_${vday}.vsdb1 $COMVSDB/${model_p}/${cmaq_domain}max_${vday}.vsdb
#  fi
fi

echo $regions

for field in $regions
do
  case $field in
    ozonemax) ver_grid="148"
              filnam="aqm.t12z.awpozcon"
              filbufr="prepbufr.tm00"
              avg_hr=24
              COMBUFR_IN1=$COMBUFR/hourly.$vday
              COMBUFR_IN2=$COMBUFR/hourly.$vdayp1
              network="ozone"
              maxgrib="maxgrib"
              XC="_ozone";;
   pmave)     ver_grid="148"
              filnam="aqm.t${cyc}z.25pm"
              filbufr="anowpm.pb.tm024"
              avg_hr=30
              COMBUFR_IN1=$COMBUFR/hourly.$vdayp1
              COMBUFR_IN2=$COMBUFR/hourly.$vdayp2
              network="pm"
              maxgrib="avegrib"
              XC="_pm";;
  pmmax)      ver_grid="148"
              filnam="aqm.t${cyc}z.25pm"
              filbufr="anowpm.pb.tm024"
              avg_hr=30
              COMBUFR_IN1=$COMBUFR/hourly.$vdayp1
              COMBUFR_IN2=$COMBUFR/hourly.$vdayp2
              network="pm"
              maxgrib="maxgrib"
              XC="_pm";;
  esac

  MODNAM=`echo $cmaq_domain | tr '[a-z]' '[A-Z]'`

  # ---------------------------------
  # Construct the keeplist parm card:
  # ---------------------------------
  rm -rf gridtobs.keeplist.${field}
  cp $PARMverf_gridtobs/verf_gridtobs.keeplist gridtobs.keeplist

  cat gridtobs.keeplist |while read line
  do
    net_name=`echo $line |awk -F"|" '{print $1}'`
    if [ $net_name = $network ]
    then
      GRID_NUM=`echo $line |awk -F"|" '{print $2}'`
      TIME_WIN=`echo $line |awk -F"|" '{print $3}'`
      OBTYP=`echo $line |awk -F"|" '{print $4}'`
      break
    fi
  done

  TIME_WIN=${vdate}

  echo "IRETGRID     - GRID NUMBER OF THE RETENTION AREA" >gridtobs.keeplist.${field}
  echo $GRID_NUM >>gridtobs.keeplist.${field}
  echo "YYMMDD       - DATE OR TIME WINDOW INDICATOR" >>gridtobs.keeplist.${field}
  echo $TIME_WIN >>gridtobs.keeplist.${field}
  echo "OBTYP        - UP TO 20 OB TYPES TO BE RETAINED" >>gridtobs.keeplist.${field}
  for obstyp in $OBTYP
  do
    echo $obstyp >>gridtobs.keeplist.${field}
  done

  # ------------------------
  # Create LEVCAT parm file:
  # ------------------------
  cp $PARMverf_gridtobs/verf_gridtobs.levcat gridtobs.levcat

  cat gridtobs.levcat |while read tmp
  do
    levcat_domain=`echo $tmp |awk -F"|" '{print $1}'`
    if [ $levcat_domain = $field ]
    then
      NUMLEV=`echo $tmp |awk -F"|" '{print $2}'`
      FIT=`echo $tmp |awk -F"|" '{print $3}'`
      break
    fi
  done

cat <<EOF_LEVCAT >gridtobs.levcat.${network}
 &LEVCAT
 NUMLEV=$NUMLEV
 FIT=$FIT
 /
EOF_LEVCAT

  typeset -Z2 fhr

  # -------------------------------
  # Start the calculation from here
  # -------------------------------
  for cmaqday in DAY1 DAY2 DAY3
  do
    case $cmaqday in
      DAY1)if [ $field = "ozonemax" ]
           then
            if [ $cyc -eq 12 ]
            then
             cp ${INDIR}.$vday/aqm.t12z.max_1hr_o3${bctag}.148.grib2 AWIP3D12.tm00_1hr_${cmaqday}
             cp ${INDIR}.$vday/aqm.t12z.max_8hr_o3${bctag}.148.grib2 AWIP3D12.tm00_8hr_${cmaqday}
             FHR=24
            fi
           elif [ $field = "pmave" -o $field = "pmmax" ]
           then
            if [ $cyc -eq 06 ]
            then
             FHR=30
             if [ $field = "pmave" ]
             then
              cp ${INDIR}.$vday/aqm.t06z.ave_24hr_pm25${bctag}.148.grib2 AWIP3D06.tm00_${field}_${cmaqday}            
             elif [ $field = "pmmax" ]
             then
              cp ${INDIR}.$vday/aqm.t06z.max_1hr_pm25${bctag}.148.grib2 AWIP3D06.tm00_${field}_${cmaqday}
             fi
            elif [ $cyc -eq 12 ]
            then
             FHR=24
             if [ $field = "pmave" ]
             then
              cp ${INDIR}.$vday/aqm.t12z.ave_24hr_pm25${bctag}.148.grib2 AWIP3D12.tm00_${field}_${cmaqday}
             elif [ $field = "pmmax" ]
             then
              cp ${INDIR}.$vday/aqm.t12z.max_1hr_pm25${bctag}.148.grib2 AWIP3D12.tm00_${field}_${cmaqday}
             fi
            fi
           fi
           ;;
    DAY2)if [ $field = "ozonemax" ]
         then
          if [ $cyc -eq 12 ]
          then
           cp ${INDIR}.$vdaym1/aqm.t12z.max_1hr_o3${bctag}.148.grib2 AWIP3D12.tm00_1hr_${cmaqday}
           cp ${INDIR}.$vdaym1/aqm.t12z.max_8hr_o3${bctag}.148.grib2 AWIP3D12.tm00_8hr_${cmaqday}
           FHR=48
          fi 
         elif [ $field = "pmave" -o $field = "pmmax" ]
         then
            if [ $cyc -eq 06 ]
            then
             FHR=54
             if [ $field = "pmave" ]
             then
              cp ${INDIR}.$vdaym1/aqm.t06z.ave_24hr_pm25${bctag}.148.grib2 AWIP3D06.tm00_${field}_${cmaqday}
             elif [ $field = "pmmax" ]
             then
              cp ${INDIR}.$vdaym1/aqm.t06z.max_1hr_pm25${bctag}.148.grib2 AWIP3D06.tm00_${field}_${cmaqday}
             fi
            elif [ $cyc -eq 12 ]
            then
             FHR=48
             if [ $field = "pmave" ]
             then
              cp ${INDIR}.$vdaym1/aqm.t12z.ave_24hr_pm25${bctag}.148.grib2 AWIP3D12.tm00_${field}_${cmaqday}
             elif [ $field = "pmmax" ]
             then
              cp ${INDIR}.$vdaym1/aqm.t12z.max_1hr_pm25${bctag}.148.grib2 AWIP3D12.tm00_${field}_${cmaqday}
             fi
          fi
         fi
         ;;
    DAY3)if [ $field = "ozonemax" ]
         then
          if [ $cyc -eq 12 ]
          then
           cp ${INDIR}.${vdaym2}/aqm.t12z.max_1hr_o3${bctag}.148.grib2 AWIP3D12.tm00_1hr_${cmaqday}
           cp ${INDIR}.${vdaym2}/aqm.t12z.max_8hr_o3${bctag}.148.grib2 AWIP3D12.tm00_8hr_${cmaqday}
           FHR=72
          fi
         elif [ $field = "pmave" -o $field = "pmmax" ]
         then
            if [ $cyc -eq 06 ]
            then
             FHR=78
             if [ $field = "pmave" ]
             then
              cp ${INDIR}.${vdaym2}/aqm.t06z.ave_24hr_pm25${bctag}.148.grib2 AWIP3D06.tm00_${field}_${cmaqday}
             elif [ $field = "pmmax" ]
             then
              cp ${INDIR}.${vdaym2}/aqm.t06z.max_1hr_pm25${bctag}.148.grib2 AWIP3D06.tm00_${field}_${cmaqday}
             fi
            elif [ $cyc -eq 12 ]
            then
             FHR=72
             if [ $field = "pmave" ]
             then
              cp ${INDIR}.${vdaym2}/aqm.t12z.ave_24hr_pm25${bctag}.148.grib2 AWIP3D12.tm00_${field}_${cmaqday}
             elif [ $field = "pmmax" ]
             then
              cp ${INDIR}.${vdaym2}/aqm.t12z.max_1hr_pm25${bctag}.148.grib2 AWIP3D12.tm00_${field}_${cmaqday}
             fi
          fi
         fi
         ;;
    esac

###    rm -rf AWIP* input_*
    rm -f input_*

   if [ $field = "ozonemax" ]
   then

    if [ $cmaqday = "DAY1" ]
     then

     if [ $cyc -eq 12 ]
     then
     $WGRIB2 AWIP3D12.tm00_1hr_${cmaqday} | grep "2147483641" | $WGRIB2 -i AWIP3D12.tm00_1hr_${cmaqday} -grib oz1hr_$cmaqday
     $WGRIB2 AWIP3D12.tm00_8hr_${cmaqday} | grep "2147483647" | $WGRIB2 -i AWIP3D12.tm00_8hr_${cmaqday} -grib oz8hr_$cmaqday
     fi
     fi

     if [ $cmaqday = "DAY2" ]
     then

     if [ $cyc -eq 12 ]
     then
     $WGRIB2 AWIP3D12.tm00_1hr_${cmaqday} | grep "17-40 hour ave" | $WGRIB2 -i AWIP3D12.tm00_1hr_${cmaqday} -grib oz1hr_$cmaqday
     $WGRIB2 AWIP3D12.tm00_8hr_${cmaqday} | grep "23-46 hour ave" | $WGRIB2 -i AWIP3D12.tm00_8hr_${cmaqday} -grib oz8hr_$cmaqday
     fi
     fi

     if [ $cmaqday = "DAY3" ]
     then

     if [ $cyc -eq 12 ]
     then
     $WGRIB2 AWIP3D12.tm00_1hr_${cmaqday} | grep "41-64 hour ave" | $WGRIB2 -i AWIP3D12.tm00_1hr_${cmaqday} -grib oz1hr_$cmaqday
     $WGRIB2 AWIP3D12.tm00_8hr_${cmaqday} | grep "47-70 hour ave" | $WGRIB2 -i AWIP3D12.tm00_8hr_${cmaqday} -grib oz8hr_$cmaqday
     fi
     fi

     cat oz1hr_$cmaqday oz8hr_$cmaqday > AWIP${avg_hr}.tm00
     ${GRB2INDEX} AWIP${avg_hr}.tm00 AWIP${avg_hr}i.tm00

cat<<eof>input_hour_max
${FHR} aqmmax     AWIP${avg_hr}.tm00
${FHR} aqmmax     AWIP${avg_hr}i.tm00
eof

cp input_hour_max input_hour

   else

      if [ $cmaqday = "DAY1" ]
      then

         if [ $cyc -eq 12 ]
         then
           $WGRIB2 AWIP3D${cyc}.tm00_${field}_${cmaqday} | grep "2147483641" | $WGRIB2 -i AWIP3D${cyc}.tm00_${field}_${cmaqday} -grib gribtemp
         else
           $WGRIB2 AWIP3D${cyc}.tm00_${field}_${cmaqday} | grep "2147483647" | $WGRIB2 -i AWIP3D${cyc}.tm00_${field}_${cmaqday}  -grib gribtemp
          fi

#      cp gribtemp AWIP3D${cyc}.tm00_${field}_${cmaqday}
       cp gribtemp ${field}_${cmaqday}
      fi

      if [ $cmaqday = "DAY2" ]
      then

         if [ $cyc -eq 12 ]
         then
           $WGRIB2 AWIP3D${cyc}.tm00_${field}_${cmaqday} | grep "17-40" | $WGRIB2 -i AWIP3D${cyc}.tm00_${field}_${cmaqday} -grib gribtemp
         else
           $WGRIB2 AWIP3D${cyc}.tm00_${field}_${cmaqday} | grep "23-46" | $WGRIB2 -i AWIP3D${cyc}.tm00_${field}_${cmaqday} -grib gribtemp
          fi

#      cp gribtemp AWIP3D${cyc}.tm00_${field}_${cmaqday}
       cp gribtemp ${field}_${cmaqday}
      fi

      if [ $cmaqday = "DAY3" ]
      then

         if [ $cyc -eq 12 ]
         then
           $WGRIB2 AWIP3D${cyc}.tm00_${field}_${cmaqday} | grep "41-64" | $WGRIB2 -i AWIP3D${cyc}.tm00_${field}_${cmaqday}  -grib gribtemp
         else
           $WGRIB2 AWIP3D${cyc}.tm00_${field}_${cmaqday} | grep "47-70" | $WGRIB2 -i AWIP3D${cyc}.tm00_${field}_${cmaqday}  -grib gribtemp
          fi

#      cp gribtemp AWIP3D${cyc}.tm00_${field}_${cmaqday}
       cp gribtemp ${field}_${cmaqday}
      fi

      cp ${field}_${cmaqday} AWIP3D${cyc}.tm00
      ${GRB2INDEX} AWIP3D${cyc}.tm00 AWIP3D${cyc}i.tm00 

cat<<eof>>input_hour
${FHR} ${field}     AWIP3D${cyc}.tm00
${FHR} ${field}     AWIP3D${cyc}i.tm00
eof

   fi

    # -----------------------------------------------
    # Create the input parm card to gridtobs program
    # -----------------------------------------------
    cat $PARMverf_gridtobs/verf_gridtobs.${field}\
    | sed -e s/MODEL_HOUR/$FHR/g -e "s/MODNAM/$MODNAM/g" -e "s/GRDNUM/$ver_grid/g" > gridtobs_${field}1
    mv gridtobs_${field}1 gridtobs_${field}_${cmaqday}

    # define a prepbufr file to filter and fit
    cp $COMBUFR_IN1/aqm.t12z.${filbufr} prepda.${vdate}
    cp $COMBUFR_IN2/aqm.t12z.${filbufr} prepda.${vdatep1}

    # ---------------------------------------------------------------
    #  run editbufr and prepfits on the combined set of observations
    # ---------------------------------------------------------------
    rm prepfits.${field}.${vdate}_$cmaqday
 
    export pgm=verf_gridtobs_editbufr
    . prep_step

#    export XLFUNIT_20=prepda.${vdate}
#    export XLFUNIT_21=prepda.${vdatep1}
#    export XLFUNIT_22=$PARMverf_gridtobs/verf_gridtobs.prepfits.tab_${network}
#    export XLFUNIT_70=bufrout_${field}

     rm -f fort.*

     ln -sf prepda.${vdate} fort.20
     ln -sf prepda.${vdatep1} fort.21
     ln -sf $PARMverf_gridtobs/verf_gridtobs.prepfits.tab_${network} fort.22
     ln -sf bufrout_${field} fort.70

    startmsg
#    $EXECverf_gridtobs/verf_gridtobs_editbufr${XC}max < gridtobs.keeplist.${field} >>editbufr_${field}_${cmaqday}.out
    $EXECverf_gridtobs/verf_gridtobs_editbufr_${field} < gridtobs.keeplist.${field} >>editbufr_${field}_${cmaqday}.out
    export err=$?; err_chk

    # ---------------------
    # Run prepfits step
    # ---------------------

    export pgm=verf_gridtobs_prepfits
    . prep_step

#    export XLFUNIT_11=gridtobs.levcat.${network}
#    export XLFUNIT_20=bufrout_${field}
#    export XLFUNIT_22=$PARMverf_gridtobs/verf_gridtobs.prepfits.tab_${network}
#    export XLFUNIT_50=prepfits.${field}.${vdate}_${cmaqday}
#    export XLFUNIT_101=airnow1_list
#    export XLFUNIT_108=airnow8_list

    rm -f fort.*

    ln -sf gridtobs.levcat.${network} fort.11
    ln -sf bufrout_${field} fort.20
    ln -sf $PARMverf_gridtobs/verf_gridtobs.prepfits.tab_${network} fort.22
    ln -sf prepfits.${field}.${vdate}_${cmaqday}${bctag} fort.50
    ln -sf airnow1_list fort.101
    ln -sf airnow8_list fort.108

    startmsg
    $EXECverf_gridtobs/verf_gridtobs_prepfits${XC} < input_hour >>prepfits_${field}_${cmaqday}${bctag}.out
    export err=$?; err_chk

    if [ $SENDCOM = YES ]
    then
      cp prepfits.${field}.${vdate}_${cmaqday} $COMOUT/aqmmax_prefits.${field}.${cmaqday}
      cp airnow1_list       ${COMOUT}/${cmaq_domain}_airnow1_list${field}.${cmaqday}
      cp airnow8_list       ${COMOUT}/${cmaq_domain}_airnow8_list${field}.${cmaqday}
      cp outgribmax ${COMOUT}/${cmaq_domain}_${field}.${cmaqday}
    fi
  
    # ----------------------
    # Run gridtobs step
    # ----------------------
    pgm=verf_gridtobs_gridtobs
    . prep_step

#    export XLFUNIT_10=prepfits.${field}.${vdate}_${cmaqday}
#    export XLFUNIT_20=$PARMverf_gridtobs/verf_gridtobs.grid104
#    export XLFUNIT_21=$PARMverf_gridtobs/verf_gridtobs.regions
#    export XLFUNIT_50=${field}_${vdate}_${cmaqday}.vdb

    rm -f fort.*

    ln -sf prepfits.${field}.${vdate}_${cmaqday}${bctag} fort.10
    ln -sf $PARMverf_gridtobs/verf_gridtobs.grid104 fort.20
    ln -sf $PARMverf_gridtobs/verf_gridtobs.regions fort.21
    ln -sf ${field}_${vdate}_${cmaqday}${bctag}.vdb fort.50

    startmsg
    $EXECverf_gridtobs/verf_gridtobs_gridtobs${XC} <gridtobs_${field}_${cmaqday} >gto_${field}_${cmaqday}${bctag}.out
    export err=$?; err_chk

    cat ${field}_${vdate}_${cmaqday}${bctag}.vdb >>${field}_${vdate}.vsdb

  done  # cmaqday

done  #field

# -------------------------------
# Save the output and vsdb files
# -------------------------------
if [ $SENDCOM = YES ]
then
  for field in $regions
  do
    cat ${field}_${vdate}.vsdb >>$COMVSDB/${model_p}/${cmaq_domain}max_${vday}.vsdb
  done
fi

#####################################################################

msg='JOB $job HAS COMPLETED NORMALLY.'
postmsg "$jlogfile" "$msg"

############## END OF SCRIPT #######################
