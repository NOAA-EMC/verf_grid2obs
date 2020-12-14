SHELL=          /bin/sh
LIBS=           $(G2_LIB4) $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
FC=             ifort
FFLAGS=         -heap-arrays -O0 
LDFLAGS=        
SRCS=        gtgdef.f  raddate.f  w3fi71.f  editbufr.f
OBJS=        gtgdef.o  raddate.o  w3fi71.o  editbufr.o
CMD=            verf_gridtobs_editbufr

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
