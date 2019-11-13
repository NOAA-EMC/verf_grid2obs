OBJS    = gtgdef.o raddate.o w3fi63.o w3fi71.o editbufr.o
SHELL=          /bin/sh
LIBS=           -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_4_64 -lbacio_4
FC=             ifort
FFLAGS=         -heap-arrays
CMD     = verf_gridtobs_editbufr_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS) 
