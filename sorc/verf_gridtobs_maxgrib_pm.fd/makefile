SHELL=          /bin/sh
LIBS=            $(W3EMC_LIB4) $(W3NCO_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
FC=             ftn
FFLAGS=         -heap-arrays -i4
OBJS	= getbak.o raddate.o prepfits.o 
CMD	= ../../exec/verf_gridtobs_maxgrib_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	rm -f *.o $(CMD)
