#!/bin/ksh
#############################################################################
# Name of script:    verf_gridtobs_fits.sh
# Purpose of script: This script generate the 24h surface and upper air 
#                    verification for various models.
# Arguments:
#  1. 'model' : Model name
#  2. 'vday' : Verification date
#  3: 'vcyc'  : model cycles
#  4. 'vlen'  : Verification length (optional)
#  5: 'frange': range of model forecast (optional) 
#############################################################################
set -x

cd $DATA

if [ $# -lt 4 ]
then
   echo "Invalid argument"
   echo "usage: verf_gridtobs.sh $model $vday $vcyc $mod"
   err_exit
fi

export model=$1
export vday=$2
export vcyc=$3
export mod=$4

###echo $mod

##export network ="GRDNAM MODNAM" 
export vdate=${vday}${vcyc}

export wgrib=${wgrib:-/gpfs/dell1/nco/ops/nwprod/grib_util.v1.1.1/exec/wgrib}
export copygb=${copygb:-/gpfs/dell1/nco/ops/nwprod/grib_util.v1.1.1/exec/copygb}

mkdir -p $DATA/$vcyc
cd $DATA/$vcyc
mkdir $DATA/$vcyc/$mod
cd $DATA/$vcyc/$mod

echo $SHELL
###typeset -L8 pll3
$utilscript/setup.sh
                                                                                         
echo $regions
for domain in $regions
do
if [ $domain = $mod ]
then
   case $domain in
     ozone) XC="_$domain"
            XCE="_$domain"
            XCT="_$domain"
            ;;
     ozonepara) XC="_ozone"
            XCE="_ozone"
            XCT="_ozone"
            ;;
     ozonepara1) XC="_ozone"
            XCE="_ozone"
            XCT="_ozone"
            ;;
     ozoneak) XC="_ozone"
            XCE="_ozone"
            XCT="_ozone"
            ;;
     ozonehi) XC="_ozone"
            XCE="_ozone"
            XCT="_ozone"
            ;;
     pm)    XC="_$domain"
            XCE="_$domain"
            XCT="_$domain"
            ;;
     pmak)  XC="_pm"
            XCE="_pm"
            XCT="_pm"
            ;;
     pmhi)  XC="_pm"
            XCE="_pm"
            XCT="_pm"
            ;;
     pm1)   XC="_pm"
            XCE="_pm"
            XCT="_pm"
            ;;
     ngac)  XC="_aod"
            XCE=""
            XCT=""
            ;;
     *)     XC=""
            XCE=""
            XCT=""
            ;;
     esac

  export DIRIN=$INDIR

  ## Get the "cycles" for the specific domain
  cp $PARMverf_gridtobs/verf_gridtobs.cycles gridtobs.cycles

  set +x
  cat gridtobs.cycles |while read line
  do
    region_name=`echo $line |awk -F"|" '{print $1}'`

    # skip the comment lines
    first_char=`echo $line |cut -c1`
    if [ "$first_char" = "#" ]
    then
       echo "It's a comment line, skip this line"
    else
      if [ $region_name = $domain ]
      then
        export model_cycles=`echo $line |awk -F"|" '{print $2}'`
        echo $model_cycles >${domain}.cycles
        break
      fi
    fi
  done
  set -x

  ## Obtain the list of model input files
  cp $PARMverf_gridtobs/verf_gridtobs.domains gridtobs.domains

#  set +x
  cat gridtobs.domains |while read line
  do
    region_name=`echo $line |awk -F"|" '{print $1}'`

    # skip the comment lines
    first_char=`echo $line |cut -c1`
    if [ "$first_char" = "#" ]
    then
       echo "It's a comment line, skip this line"
    else
      if [ $region_name = $domain ]
      then
        export network=`echo $line |awk -F"|" '{print $2}'`
        export runnam=`echo $line |awk -F"|" '{print $3}'`
        export MODNAM=`echo $line |awk -F"|" '{print $4}'`
        export pll3=`echo $line |awk -F"|" '{print $5}'`
        export sfhr=`echo $line |awk -F"|" '{print $6}'`
        export clen=`echo $line |awk -F"|" '{print $7}'`
        export inc=`echo $line |awk -F"|" '{print $8}'`
        export filnam1=`echo $line |awk -F"|" '{print $9}'`
        export filnam2=`echo $line |awk -F"|" '{print $10}'`
        export filnam3=`echo $line |awk -F"|" '{print $11}'`
        export tmkk=`echo $line |awk -F"|" '{print $12}'`
        export DOCYGB=`echo $line |awk -F"|" '{print $13}'`
        export cgrid=`echo $line |awk -F"|" '{print $14}'`
        export GRDNAM=`echo $line |awk -F"|" '{print $15}'`
        export GRIB=`echo $line |awk -F"|" '{print $16}'`

        break
      fi
    fi
  done
#  set -x

  $USHverf_gridtobs/verf_gridtobs_getmodinput.sh $domain

  ##  Get the prepbufr files 
#  if [ $model = namx ]
#  then
#  $USHverf_gridtobs/verf_gridtobs_getbufr_retro.sh $vdate $domain
#  else
  . $USHverf_gridtobs/verf_gridtobs_getbufr.sh $vdate $domain
#  fi

##  echo "here??"
#  chgrp rstprod prepda.${vdate} 
  
##  . /u/Jack.Woollen/bin/gbq prepda.${vdate}
##    /nwprod/util/exec/cwordsh unblk prepda.${vdate} prepda.new

    chgrp rstprod prepda.${vdate}

#   echo "before keeplist"
#   echo $network

  ## Build the keeplist parm card from template:
  if [ ! -s gridtobs.keeplist.${network} ]
  then
    cp $PARMverf_gridtobs/verf_gridtobs.keeplist gridtobs.keeplist

    set +x
    cat gridtobs.keeplist |while read line
    do
      net_name=`echo $line |awk -F"|" '{print $1}'`
      if [ $net_name = $network ]
      then
        GRID_NUM=`echo $line |awk -F"|" '{print $2}'` 
        if [ $GRID_NUM = "0" ] 
        then
         LAT_LON=`echo $line |awk -F"|" '{print $3}'`
         TIME_WIN=`echo $line |awk -F"|" '{print $4}'`
         OBTYP=`echo $line |awk -F"|" '{print $5}'`
        else
         TIME_WIN=`echo $line |awk -F"|" '{print $3}'`
         OBTYP=`echo $line |awk -F"|" '{print $4}'`
        fi
    
        break
      fi

    done

    ## Special treatment for the AQM model:
    if [ $model = "aqm" -o $model = "aqmak" -o $model = "pm" -o $model = "aqmpara1" -o $model = "pmak" -o $model = "pm1" -o $model = "aqmhi" -o $model = "pmhi" -o $model = "aqmpara" ]
    then
       vdatep1=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.2/exec/ips/ndate +24 $vdate`
       if [ $vcyc = 24 ]
       then 
          TIME_WIN=${vdatep1}
       else
          TIME_WIN=${vdate}
       fi
    fi 
 
    set -x

    rm -rf gridtobs.keeplist.${network}
    
    echo "IRETGRID     - GRID NUMBER OF THE RETENTION AREA" >gridtobs.keeplist.${network}
    echo $GRID_NUM >>gridtobs.keeplist.${network}
    if [ $GRID_NUM = "0" ] 
    then
     echo $LAT_LON >>gridtobs.keeplist.${network}
    fi
    echo "YYMMDD       - DATE OR TIME WINDOW INDICATOR" >>gridtobs.keeplist.${network}
    echo $TIME_WIN >>gridtobs.keeplist.${network}
    echo "OBTYP        - UP TO 20 OB TYPES TO BE RETAINED" >>gridtobs.keeplist.${network}
    for obstyp in $OBTYP
    do
       echo $obstyp >>gridtobs.keeplist.${network}
    done
  fi

  export pgm=verf_gridtobs_editbufr${XCE}
  . $utilscript/prep_step

   pwd
   echo $DATA

#  export XLFUNIT_20=prepda.${vdate}
#  export XLFUNIT_50=bufrout_${domain}

   ln -sf $DATA/$vcyc/$mod/prepda.${vdate} fort.20
   ln -sf $DATA/$vcyc/$mod/bufrout_${domain} fort.50

  $utilscript/startmsg.sh
  $EXECverf_gridtobs/verf_gridtobs_editbufr${XCE} < gridtobs.keeplist.${network} >>$pgmout
  export err=$?; $utilscript/err_chk.sh

  rm fort.*

  ## Create LEVCAT parm file:
  cp $PARMverf_gridtobs/verf_gridtobs.levcat gridtobs.levcat

  set +x
  cat gridtobs.levcat |while read tmp
  do
    levcat_domain=`echo $tmp |awk -F"|" '{print $1}'`
    if [ $levcat_domain = $domain ]
    then
       NUMLEV=`echo $tmp |awk -F"|" '{print $2}'`
       FIT=`echo $tmp |awk -F"|" '{print $3}'`
       break
    fi
  done
#  set -x

  cat <<EOF_LEVCAT >gridtobs.levcat.${domain}
 &LEVCAT
 NUMLEV=$NUMLEV
 FIT=$FIT
 /
EOF_LEVCAT

  ## Run prepfits to comine the bufr data and the model grid output data
  for datacard in `ls input_${domain}*`
  do
    filesize=`cat nonzero_${domain}|wc -l`
    if [ $filesize -le 0 ]; then
       echo "Warning: Input file is empty, no verification for model $domain"
    else
       export pgm=verf_gridtobs_prepfits
       . $utilscript/prep_step

#       export XLFUNIT_11=gridtobs.levcat.${domain}
#       export XLFUNIT_20=bufrout_${domain}
#       export XLFUNIT_22=$PARMverf_gridtobs/verf_gridtobs.prepfits.tab${XCT}
#       export XLFUNIT_50=prepfits.${domain}.${vdate}

       ln -sf $DATA/$vcyc/$mod/gridtobs.levcat.${domain} fort.11
       ln -sf $DATA/$vcyc/$mod/bufrout_${domain} fort.20
       ln -sf $PARMverf_gridtobs/verf_gridtobs.prepfits.tab${XCT} fort.22
       ln -sf $DATA/$vcyc/$mod/prepfits.${domain}.${vdate} fort.50

       $utilscript/startmsg.sh
       echo "datacard"
       echo $datacard
       $EXECverf_gridtobs/verf_gridtobs_prepfits${XC} < $datacard >>prepfit.out.${domain}
       export err=$?
       cat prepfit.out.${domain} >>$pgmout
       $utilscript/err_chk.sh

       ## Save the "prepfits" files 
       chmod 640 prepfits.${domain}.${vdate}
       chgrp rstprod prepfits.${domain}.${vdate}
       cp prepfits.${domain}.${vdate} $COMOUT/${model}.t${vcyc}z.prepfits.${domain}
#      chgrp rstprod $COMOUT/${model}.t${vcyc}z.prepfits.${domain}
       if [ $domain = "srmean" -o $domain = "srmeanx" ]
       then
         cp hour_${domain} $COMOUT/hour_${domain}.${vdate}
       fi

       if [ $domain = "gec00" ]
       then
         cp hour_${domain} $COMOUT/hour_${domain}.${vdate}
       fi
       
       ## Now run gridtobs to generate the statistics
       rm -rf gridtobs_${domain}

       ## Copy in the generic parameter file:
       cp $PARMverf_gridtobs/verf_gridtobs.${network}  verf_gridtobs.${domain}

#       if [ $model = "ruc" ]
#        then
#          cp $PARMverf_gridtobs/verf_gridtobs.${network}2 verf_gridtobs.${domain}
#       fi

       ## Create the input parm card by inserting the hours:
       nhours=`cat hour_${domain} |wc -l`
       shr=`sed -n 1p hour_${domain}`

       nline=`cat verf_gridtobs.${domain} |wc -l`
       i=1
       jstat=0

       echo $SHELL
       typeset -R5 nhrs
       nhrs=$nhours
  
       set +x
       while [ $i -le $nline ]
       do
         jsr=`sed -n ${i}p verf_gridtobs.${domain}  |awk '{print $1}'`
         if [ "$jsr" = "V01" ]; then  
            let "jstat=jstat+1"
         fi
 
        ## Skip the first V01 section in the parameter file for the off-time cycles
         if [[ `expr $vcyc % 12` -eq 0 ]] || [[ $jstat -gt 1 ]] || [[ $model = "aqm" ]] || [[ $model = "aqmak" ]] || [[ $model = "rtma" ]] || [[ $model = "satimg" ]] || [[ $model = "satimgb" ]] || [[ $model = "satimgx" ]]|| [[ $model = "ngac" ]] || [[ $model = "pm" ]] || [[ $model = "aqmpara1" ]] ||  [[ $model = "pm1" ]] ||  [[ $model = "pmak" ]] ||  [[ $model = "aqmhi" ]] ||  [[ $model = "pmhi" ]] || [[ $model = "aqmpara" ]]
         then
           sed -n ${i}p verf_gridtobs.${domain} >>gridtobs_${domain}
           isr=`sed -n ${i}p verf_gridtobs.${domain}  |awk -F"/" '{print $2}'`
           if [ "$isr" = "GRDNAM" ]
           then
             echo "$nhrs  $shr" >>gridtobs_${domain}

             # Comment out this section--20080731
             #let "hh=shr+inc"
             #while [ $hh -le $frange ]
             #do
             #  if [ $hh -lt 10 ]; then hh=0$hh; fi
             #  echo "  $hh">>gridtobs_${domain}
             #  let "hh=hh+inc"
             #done
             cat hour_${domain} |while read hh
             do
               if [ $hh -eq $shr ]
               then
                 echo "Skip the first hour"
               else
                 echo "  $hh">>gridtobs_${domain}
               fi
             done
            
           fi
         fi
         let "i=i+1"
       done  
       set -x 

       sed -e "s/MODNAM/$MODNAM/g" -e "s/GRDNAM/$GRDNAM/g" gridtobs_${domain} >gridtobs_${domain}.new
       mv gridtobs_${domain}.new gridtobs_${domain}

#       cat firemask00.grb firemask12.grb firemask24.grb firemask36.grb > firemask.grb
#       $utilexec/grbindex firemask.grb firemaski.grb       
 
       pgm=verf_gridtobs_gridtobs
       . $utilscript/prep_step

#       export XLFUNIT_10=prepfits.${domain}.${vdate}
#       export XLFUNIT_20=$PARMverf_gridtobs/verf_gridtobs.grid104
#       export XLFUNIT_21=$PARMverf_gridtobs/verf_gridtobs.regions
#       export XLFUNIT_50=${domain}_${vdate}.vdb

       rm -f fort.*

       ln -sf $DATA/$vcyc/$mod/prepfits.${domain}.${vdate} fort.10
       ln -sf $PARMverf_gridtobs/verf_gridtobs.grid104 fort.20
       ln -sf $PARMverf_gridtobs/verf_gridtobs.regions fort.21
       ln -sf $DATA/$vcyc/$mod/${domain}_${vdate}.vdb fort.50

       $utilscript/startmsg.sh
       $EXECverf_gridtobs/verf_gridtobs_gridtobs${XC} <gridtobs_${domain} >gto.${domain}${vcyc}.out
       export err=$?; $utilscript/err_chk.sh

       cat gto.${domain}${vcyc}.out >>$pgmout

       cat ${domain}_${vdate}.vdb >>${domain}_${vday}${vcyc}.vsdb
#
# REPEAT GRIDTOBS FOR VERIF ON FIRE WX DOMAINS FOR NAM AND RUC - 5 Apr 2011 - PS
#
       if [ $model = "nam" -o $model = "namak" -o $model = "rap" -o $model = "rapak" -o $model = "conusnest" -o $model = "aknest" -o $model = "hrrr" -o $model = "hiresw" -o $model = "nssl4arw" ]
         then
          cp $PARMverf_gridtobs/verf_gridtobs.${network}_fwis  verf_gridtobs.${domain}_fwis 

## Create the input parm card by inserting the hours:
       nhours=`cat hour_${domain}_fwis |wc -l`
       shr=`sed -n 1p hour_${domain}_fwis`

       nline=`cat verf_gridtobs.${domain}_fwis |wc -l`
       i=1
       jstat=0

       typeset -R5 nhrs
       nhrs=$nhours
  
       set +x
       while [ $i -le $nline ]
       do
         jsr=`sed -n ${i}p verf_gridtobs.${domain}_fwis  |awk '{print $1}'`
         if [ "$jsr" = "V01" ]; then  
            let "jstat=jstat+1"
         fi

#         if [[ `expr $vcyc % 12` -eq 0 ]] || [[ $jstat -gt 1 ]]
#           then
           sed -n ${i}p verf_gridtobs.${domain}_fwis >>gridtobs_${domain}_fwis
           isr=`sed -n ${i}p verf_gridtobs.${domain}  |awk -F"/" '{print $2}'`
           if [ "$isr" = "GRDNAM" ]
           then
             echo "$nhrs  $shr" >>gridtobs_${domain}_fwis

             cat hour_${domain}_fwis |while read hh
             do
               if [ $hh -eq $shr ]
               then
                 echo "Skip the first hour"
               else
                 if [ $hh -le 36 ]
                 then
                   echo "  $hh">>gridtobs_${domain}_fwis
                 fi
               fi
             done
            
           fi
#         fi
         let "i=i+1"
       done  
       set -x 

       sed -e "s/MODNAM/$MODNAM/g" -e "s/GRDNAM/$GRDNAM/g" gridtobs_${domain}_fwis >gridtobs_${domain}_fwis.new
       mv gridtobs_${domain}_fwis.new gridtobs_${domain}_fwis

       echo $vcyc > cycfile
       cat cycfile gridtobs_${domain}_fwis > gridtobs_${domain}_fwis.new
       mv gridtobs_${domain}_fwis.new gridtobs_${domain}_fwis

       pgm=verf_gridtobs_gridtobs_fwis
       . $utilscript/prep_step

#       export XLFUNIT_10=prepfits.${domain}.${vdate}
#       export XLFUNIT_20=$PARMverf_gridtobs/verf_gridtobs.grid104
#       export XLFUNIT_21=$PARMverf_gridtobs/verf_gridtobs.regions
#       export XLFUNIT_31=firemask1.grb
#       export XLFUNIT_32=firemask1i.grb
#       export XLFUNIT_33=firemask2.grb
#       export XLFUNIT_34=firemask2i.grb
#       export XLFUNIT_35=firemask3.grb
#       export XLFUNIT_36=firemask3i.grb
#       export XLFUNIT_37=firemask4.grb
#       export XLFUNIT_38=firemask4i.grb
#       export XLFUNIT_39=firemask5.grb
#       export XLFUNIT_40=firemask5i.grb
#       export XLFUNIT_41=firemask6.grb
#       export XLFUNIT_42=firemask6i.grb
#       export XLFUNIT_43=firemask7.grb
#       export XLFUNIT_44=firemask7i.grb
#       export XLFUNIT_50=${domain}_${vdate}_fwis.vdb

        rm fort.*
        ln -sf prepfits.${domain}.${vdate} fort.10
        ln -sf $PARMverf_gridtobs/verf_gridtobs.grid104 fort.20
        ln -sf $PARMverf_gridtobs/verf_gridtobs.regions fort.21
        ln -sf firemask1.grb fort.31
        ln -sf firemask1i.grb fort.32
        ln -sf firemask2.grb fort.33
        ln -sf firemask2i.grb fort.34
        ln -sf firemask3.grb fort.35
        ln -sf firemask3i.grb fort.36
        ln -sf firemask4.grb fort.37
        ln -sf firemask4i.grb fort.38
        ln -sf firemask5.grb fort.39
        ln -sf firemask5i.grb fort.40
        ln -sf firemask6.grb fort.41
        ln -sf firemask6i.grb fort.42
        ln -sf firemask7.grb fort.43
        ln -sf firemask7i.grb fort.44
        ln -sf ${domain}_${vdate}_fwis.vdb fort.50
       $utilscript/startmsg.sh
       $EXECverf_gridtobs/verf_gridtobs_gridtobs_fwis <gridtobs_${domain}_fwis >gto.${domain}${vcyc}_fwis.out
       export err=$?; $utilscript/err_chk.sh

       cat gto.${domain}${vcyc}_fwis.out >>$pgmout

       cat ${domain}_${vdate}_fwis.vdb >>${domain}_${vday}${vcyc}.vsdb
      fi
## End of Fire wx on NAM/RUC verification

    fi
   
    
  done

  # Copy vsdb files to /com and alert them to the server if needed:
  modnam=`echo $MODNAM | tr '[A-Z]' '[a-z]'`

  if [ $SENDCOM = YES ]
  then
    cat $DATA/$vcyc/$mod/${domain}_${vday}${vcyc}.vsdb >>$COMVSDB/${model}/${modnam}_${vday}.vsdb
  fi

fi
done

## DBN alert the vsdb record files to TOC:
if [ $SENDDBN = YES ]
then
   cd $COMVSDB/$model
   for vsdbfile in `ls $model/*_${vday}.vsdb`
   do 
     $DBNROOT/bin/dbn_alert MODEL VERIF_GRID2OBS $job $COMVSDB/${model}/$vsdbfile
   done
fi
