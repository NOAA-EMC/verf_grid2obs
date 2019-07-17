SHELL=          /bin/sh
LIBS =       -L/nwprod/lib -lg2_4 -ljasper -lw3nco_4 -lw3emc_4 -lbufr_s_64 -lbacio_4 -lpng -lz
FC=             ifort
INC =        /nwprod/lib/incmod/g2_4
FFLAGS =     -i4 -I$(INC)
LDFLAGS=       
SRCS=       writemask.f
OBJS=       writemask.o 
CMD=        verf_gridtobs_writemask

all:    $(CMD)

$(CMD): $(OBJS)
	$(FC) $(FFLAGS)$(LDFLAGS) -g  -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
