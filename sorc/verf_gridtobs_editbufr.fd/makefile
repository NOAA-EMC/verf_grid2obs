SHELL=          /bin/sh
LIBS=           -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_4_64 -lbacio_4
FC=             ifort
FFLAGS=         -heap-arrays
LDFLAGS=
SRCS=        gtgdef.f  raddate.f  w3fi71.f  editbufr.f
OBJS=        gtgdef.o  raddate.o  w3fi71.o  editbufr.o
CMD=            verf_gridtobs_editbufr

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
