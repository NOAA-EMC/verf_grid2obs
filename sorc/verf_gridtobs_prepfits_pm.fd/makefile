OBJS	= calcape.o adddate.o getbak.o getfct.o getprof.o ireadpb.o raddate.o \
	setterp.o ufbxy3.o w3fb13.o prepfits.o
SRCS	= calcape.f adddate.f getbak.f getfct.f getprof.f ireadpb.f raddate.f\
	setterp.f ufbxy3.f w3fb13.f prepfits.f
SHELL=          /bin/sh
LIBS =       -L/nwprod/lib -lg2_4 -ljasper -lw3nco_4 -lw3emc_4 -lbufr_s_64 -lbacio_4 -lpng -lz
INC =        /nwprod/lib/incmod/g2_4
FFLAGS =     -i4 -g -traceback -I$(INC)
LDFLAGS =    -shared-intel
FC=             ifort
CMD	= verf_gridtobs_prepfits_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	rm *.o
