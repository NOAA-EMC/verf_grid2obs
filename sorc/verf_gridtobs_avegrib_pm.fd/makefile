SHELL=          /bin/sh
LIBS=            $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIBs) $(BACIO_LIB4)
FC=             ifort
FFLAGS=         -heap-arrays
OBJS	= getbak.o raddate.o prepfits.o 
CMD	= verf_gridtobs_avegrib_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	rm -f *.o $(CMD)

