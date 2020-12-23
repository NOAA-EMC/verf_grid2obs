SHELL=          /bin/sh
LIBS =        $(G2_LIB4) $(JASPER_LIB) $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIB4) $(BACIO_LIB4) $(PNG_LIB) $(Z_LIB)
FC=             ifort
INC =        $(G2_INC4)
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
