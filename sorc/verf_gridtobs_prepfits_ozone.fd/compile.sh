rm prep.a
ncepxlf -qcheck -qextchk -qfixed -c *.f
/bin/rm prepfits.o
ar -ruv -X 64 prep.a *.o
ncepxlf -qcheck -qextchk -qfixed -g -o verif_ozone_max_prepfits \
        prepfits.f prep.a \
      /nwprod/lib/libw3_4.a \
      /nwprod/lib/libbufr_4_64.a \
      /nwprod/lib/libbacio_4.a -bmaxdata:2000000000 -bmaxstack:256000000
      
rm *.o
