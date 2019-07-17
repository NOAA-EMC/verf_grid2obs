SHELL=          /bin/sh
LIBS=           -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_s_64 -lbacio_4
FC=             ifort
FFLAGS=         -heap-arrays
OBJS	= getbak.o raddate.o prepfits.o 
CMD	= verf_gridtobs_avegrib_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

