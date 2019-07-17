SHELL=          /bin/sh
LIBS=           -L/nwprod/lib -lw3emc_4 -lw3nco_4 -lbufr_4_64 -lbacio_4
FC=             ifort
FFLAGS=         -heap-arrays -i4
OBJS	= getbak.o raddate.o prepfits.o 
CMD	= verf_gridtobs_maxgrib_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

