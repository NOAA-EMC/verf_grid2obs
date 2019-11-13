#!/bin/ksh

#####################################################################
# Script to run grid2obs job for SREF grid2obs verification on 
# Eta, EKF, RSM, WRF and ctl groups and output in vsdb record. 
# part I: ANYAIR, PROFLR, VADWND
#
# Fields verified (individual members):
#  Z/T/Wind (11 p-levs: 1000,850,700,500,400,300,250,200,150,100,50mb)
#  RH (6 p-levs:        1000,850,700,500,400,300mb)
#  surface:             SLP, T2m, RH2m, 10m-vector-wind
#
# Log
#  Before May 2006: Michael Baker
#  Jan/Feb, 2007:   Jun Du, combined the original 80 small scripts into
#                   one and extended from 63hr to 87hr
#  ???, 2007:       combine with Eric Rogers's package to include ensemble mean??
#  April, 2008:     Julia Zhu, Reformat the script for production implementation 
#################################################################################
 
set -x

echo $SHELL

export PARMverf_gridtobs=${PARMverf_gridtobs:-/meso/save/wx20ps/verif/nwprod/parm}
export EXECverf_gridtobs=${EXECverf_gridtobs:-/meso/save/wx20ps/verif/nwprod/exec}
export frange=${frange:-87}
export inc=${inc:-06}

if [ $# -lt 1 ]
then
   echo "Invalid argument"
   echo "usage: verf_gridtobs_ens_sref.sh $domain"
   err_exit
fi

domain=$1

case $domain in
  srefens)  export XC="";
            export MODEL=SREFENS;;
  srefensx) export XC="x";
            export MODEL=SREFENSX;;   
esac

####rm -rf fort.*
###for group in ctl all
for group in ctl arw nmm nmb ens
do
#  case $group in
#    ctl) ens_mem="ebmjc${XC} ekfcc${XC} rsasc${XC} emc${XC} nmmc${XC}"
#         integer NUMFIL=30
#         ;;
#    all) ens_mem="ebmjc${XC} ebmjp${XC} ebmjn${XC} ekfcc${XC} ekfcp${XC} ekfcn${XC} \
#                  rsasc${XC} rsasp${XC} rsasn${XC} rrasp${XC} rrasn${XC} \
#                  emc${XC} emp1${XC} emn1${XC} emp2${XC} emn2${XC} \
#                  nmmc${XC} nmmp1${XC} nmmn1${XC} nmmp2${XC} nmmn2${XC}"
#         integer NUMFIL=60
#         ;;
#  esac
#
#  for mem in $ens_mem
#  do
#     if [ -e $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} ]
#     then
#      ln -sf $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem}   fort.${NUMFIL}
#      ((NUMFIL = NUMFIL + 1))
#     fi
#  done
   if [ $domain = "srefens" ] 
   then
     if [ $vcyc -eq 0 ] || [ $vcyc -eq 12 ]
     then
       bytemin=2000000
     else  
       bytemin=500000
     fi
   elif [ $domain = "srefensx" ]
   then
     if [ $vcyc -eq 0 ] || [ $vcyc -eq 12 ]
     then
       bytemin=1500000
     else
       bytemin=400000
     fi
   fi
   case $group in
     ctl) ens_mem="srarwc${XC} srnmbc${XC}"
          integer NUMFIL=30
          integer mod_num_ctl=0
          for mem in $ens_mem
           do
             if [ -e $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} ]
             then
             a=`ls -l $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} | awk '{ print $5 }'`
             if [ $a -gt $bytemin ]
              then
               ln -sf $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem}   fort.${NUMFIL}
               MEM=`echo $mem | tr a-z A-Z`
               if [ $NUMFIL -eq 30 ] 
               then
                 first_mem_ctl=$MEM

                 ENS_MEM_CTL="$ENS_MEM_CTL $MEM"
               fi
               echo $mem $MEM
               echo $first_mem_ctl
               echo $ENS_MEM_CTL
               ((NUMFIL = NUMFIL + 1))
               ((mod_num_ctl = mod_num_ctl + 1))
             fi
             fi
           done
           MODNAM_CTL=`echo $group | tr a-z A-Z`
           echo $MODNAM_CTL
           ver_num_ctl=30
           gname_ctl="ctl"
           ;;
     arw) ens_mem="srarwc${XC} srarwp1${XC} srarwp2${XC} srarwp3${XC} srarwp4${XC} srarwp5${XC} srarwp6${XC} srarwn1${XC} srarwn2${XC} srarwn3${XC} srarwn4${XC} srarwn5${XC} srarwn6${XC} "
          integer NUMFIL=60
          ver_num_arw=$NUMFIL
          integer mod_num_arw=0
          for mem in $ens_mem
          do
             if [ -e $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} ]
             then
             a=`ls -l $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} | awk '{ print $5 }'`
             if [ $a -gt $bytemin ]
              then
               ln -sf $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem}   fort.${NUMFIL}
               MEM=`echo $mem | tr a-z A-Z`
               if [ $NUMFIL -eq 60 ]
               then
                 first_mem_arw=$MEM
               else
                 ENS_MEM_ARW="$ENS_MEM_ARW $MEM"
               fi
               echo $mem $MEM
               echo $first_mem_arw
               echo $ENS_MEM_ETA
               ((NUMFIL = NUMFIL + 1))
               ((mod_num_arw = mod_num_arw + 1))
             fi
             fi
           done
           MODNAM_ARW=`echo $group | tr a-z A-Z`
           echo $MODNAM_ARW
           NUMSAVE=$NUMFIL
           gname_arw="all"
           ;;
       nmb) ens_mem="srnmbc${XC} srnmbp1${XC} srnmbp2${XC} srnmbp3${XC} srnmbp4${XC} srnmbp5${XC} srnmbp6${XC} srnmbn1${XC} srnmbn2${XC} srnmbn3${XC} srnmbn4${XC} srnmbn5${XC} srnmbn6${XC}"
            integer NUMFIL=$NUMSAVE
            ver_num_nmb=$NUMFIL
            integer mod_num_nmb=0
            for mem in $ens_mem
          do
             if [ -e $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} ]
             then
             a=`ls -l $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem} | awk '{ print $5 }'`
             if [ $a -gt $bytemin ]
              then
               ln -sf $COMIN/sref${XC}.t${vcyc}z.prepfits.${mem}   fort.${NUMFIL}
               MEM=`echo $mem | tr a-z A-Z`
               if [ $NUMFIL -eq $NUMSAVE ]
               then
                 first_mem_nmb=$MEM
               else
                 ENS_MEM_NMB="$ENS_MEM_NMB $MEM"
               fi
               echo $mem $MEM
               echo $first_mem_nmb
               echo $ENS_MEM_NMB
               ((NUMFIL = NUMFIL + 1))
               ((mod_num_nmb = mod_num_nmb + 1))
             fi
             fi
           done
           MODNAM_NMB=`echo $group | tr a-z A-Z`
           echo $MODNAM_NMB
           NUMSAVE=$NUMFIL
           gname_nmb="all"
           ;;
        ens) first_mem_ens=$first_mem_arw
             ENS_MEM_ENS="$ENS_MEM_ARW $first_mem_nmm $ENS_MEM_NMM $first_mem_nmb $ENS_MEM_NMB"
             ((mod_num_ens = mod_num_arw + mod_num_nmm + mod_num_nmb))
             ver_num_ens=60
             MODNAM_ENS=`echo $group | tr a-z A-Z`
             gname_ens="ens"
     esac
          
done

rm -rf srefensp${XC}_${vdate}.vsdb

## Now run gridtobs_ens to generate the VSDB records
CYEAR=`echo $vday |cut -c1-4`
SUMMER_S=${CYEAR}0323
SUMMER_E=${CYEAR}0923

if [ $vday -lt $SUMMER_S -o $vday -gt $SUMMER_E ]
then
  cp $PARMverf_gridtobs/verf_gridtobs.srefens.winter  verf_gridtobs.${domain}
else
  cp $PARMverf_gridtobs/verf_gridtobs.srefens.summer  verf_gridtobs.${domain}
fi

rm -rf poe_ens

NUM_PROC=0

#################################################################
# Prepare the input parm card and generate the VSDB files 
# Note: Divide the ENS (21 members) into 2 parts according to
#       the hours due to the long wall-o'clock time it needs
#################################################################

## Create the input parm card by inserting the hours:
cp $COMIN/hour_srmean${XC}.${vdate} hour_${domain}
#shr=`sed -n 1p hour_${domain}`
#ehr=`tail -1 hour_${domain}`
nhours=`cat hour_${domain}|wc -l`

## Divide the number of hours into 2 pieces:
rm -rf hourlist_1 hourlist_2 hourlist_3
j=0
for hh in `cat hour_${domain}`
do
  if [ `expr $j % 3` -eq 0 ]
  then
    echo $hh >>hourlist_1
  elif [ `expr $j % 3` -eq 1 ]
  then
    echo $hh >>hourlist_2
  else 
    echo $hh >>hourlist_3
  fi
  let "j=j+1"
done

rm -rf gridtobs_${domain}_*

jseg=01

for list in 1 2 3
do

  ###################################
  # Construct the parm file
  ###################################

  echo $iseg $jseg $nhrs $mod_num
  echo $SHELL
  typeset -Z2 iseg jseg
  typeset -R5 nhrs mod_num
  nhrs=`cat hourlist_${list}|wc -l`
  shr=`sed -n 1p hourlist_${list}`
  ehr=`tail -1 hourlist_${list}`

  for group in ctl arw nmb ens
  do
    case $group in
     ctl) MODNAM=$MODNAM_CTL
          first_mem=$first_mem_ctl
          ens_mem=$ENS_MEM_CTL 
          GRDNAM=212
          mod_num=$mod_num_ctl
          ver_num=$ver_num_ctl
          gname=$gname_ctl
          ;;
     arw) MODNAM=$MODNAM_ARW
          first_mem=$first_mem_arw
          ens_mem=$ENS_MEM_ARW
          GRDNAM=212
          mod_num=$mod_num_arw
          ver_num=$ver_num_arw
          gname=$gname_arw
          ;;
     nmb) MODNAM=$MODNAM_NMB
          first_mem=$first_mem_nmb
          ens_mem=$ENS_MEM_NMB 
          GRDNAM=212
          mod_num=$mod_num_nmb
          ver_num=$ver_num_nmb
          gname=$gname_nmb
          ;;
     ens) MODNAM=$MODNAM_ENS
          first_mem=$first_mem_ens
          ens_mem=$ENS_MEM_ENS
          GRDNAM=212
          mod_num=$mod_num_ens
          ver_num=$ver_num_ens
          gname=$gname_ens
          ;;
    esac

    jstat=0
    iseg=0

    nline=`cat verf_gridtobs.${domain} |wc -l`
    iline=1

    while [ $iline -le $nline ]
    do
      isr=`sed -n ${iline}p verf_gridtobs.${domain}  |awk '{print $1}'`

      if [ "$isr" = "V01" ]
      then
        jstat=`expr $jstat + 1`
      fi
    
      if [[ `expr $vcyc % 12` -eq 0  ]] || [[ $jstat -gt 1 ]]
      then
        if [ "$isr" = "V01" ]
        then
          iseg=`expr $iseg + 1`

          if [ $gname = "ctl" ]
          then
            inum=$jseg
          else
            inum=$iseg
          fi
          
          sed -n ${iline}p verf_gridtobs.${domain} >>gridtobs_${domain}_${gname}_${list}_$inum
          sed -e "s/MODEL/$MODEL/g" -e "s/MODNAM/$MODNAM/g" -e "s/VER_NUM/$ver_num/g" \
                 gridtobs_${domain}_${gname}_${list}_$inum >gridtobs_${domain}_${gname}_${list}_${inum}.1
          mv gridtobs_${domain}_${gname}_${list}_${inum}.1 gridtobs_${domain}_${gname}_${list}_$inum
          echo "$mod_num  $first_mem/$GRDNAM" >>gridtobs_${domain}_${gname}_${list}_$inum
          for mem in $ens_mem
          do
            echo "  $mem/$GRDNAM" >>gridtobs_${domain}_${gname}_${list}_$inum
          done 
          echo "$nhrs  $shr" >>gridtobs_${domain}_${gname}_${list}_$inum
          for hh in `cat hourlist_${list}`
          do
            if [ $hh -ne $shr ]
            then
              echo "  $hh">>gridtobs_${domain}_${gname}_${list}_$inum
            fi
          done
        else
          sed -n ${iline}p verf_gridtobs.${domain} >>gridtobs_${domain}_${gname}_${list}_$inum
        fi
      fi

     let "iline=iline+1"
    done   # iline

  done  #group

  let "jseg=jseg+1"
done   # list

## End of generating the parm card

#ln -sf $PARMverf_gridtobs/verf_gridtobs.grid104 fort.20
#ln -sf $PARMverf_gridtobs/verf_gridtobs.regions fort.21

ulimit -s unlimited

for parmcard in `ls gridtobs_${domain}_???_*`
do
  echo "$USHverf_gridtobs/verf_gridtobs_srefens.sh $parmcard">>poe_ens
  ## Total Number of processors needed:
  NUM_PROC=`expr $NUM_PROC + 1`
done

# Execute the poe script:

mv poe_ens cmdfile
chmod +x cmdfile

#export MP_CMDFILE=cmdfile
#export MP_EUIDEVICE=sn_all
#export MP_EUILIB=us
#export MP_PGMMODEL=mpmd
#export MP_TASK_AFFINITY=cpu
#export MP_LABELIO=YES
#export MP_INFOLEVEL=3
#export MP_PULSE=0

#cp /u/Farid.Parpia/trip/testmalloc .
#mpirun.lsf /usr/bin/time -f "max RSS (times 4): %M KB" ./testmalloc 1024

#mpirun.lsf
mpirun cfp cmdfile


export err=$?; err_chk

cat ${domain}_*_${vdate}.vdb >>${domain}_${vdate}.vsdb

if [ $SENDCOM = YES ]
then
  if [ ! -s $COMVSDB/${domain} ]; then mkdir -p $COMVSDB/${domain}; fi
  if [ -s $COMVSDB/${domain}/${domain}_${vday}.vsdb ]
  then
    sed -e "/$vday$vcyc/d" $COMVSDB/${domain}/${domain}_${vday}.vsdb >$COMVSDB/${domain}/${domain}_${vday}.vsdb1
    mv $COMVSDB/${domain}/${domain}_${vday}.vsdb1 $COMVSDB/${domain}/${domain}_${vday}.vsdb
    cat ${domain}_${vdate}.vsdb >>$COMVSDB/${domain}/${domain}_${vday}.vsdb
  else
    cp ${domain}_${vdate}.vsdb $COMVSDB/${domain}/${domain}_${vday}.vsdb
  fi

fi
