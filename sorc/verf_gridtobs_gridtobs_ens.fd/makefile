SHELL=          /bin/sh
LIBS=            $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
FC=             ifort
FFLAGS=         -g -pg -O -xHost -i4
LDFLAGS=       
SRCS=       getnumber.f gtgdef.f igotdata.f inarea.f inlayer.f isitob.f readcntl.f \
            setarea.f setfcst.f setlevel.f setmodel.f setobtyp.f stclst.f stldsp.f \
            stlstr.f stnumb.f stread.f strmbl.f strxbl.f \
            w3fb06.f w3fb11.f w3fb71.f gridtobs.f
OBJS=       getnumber.o gtgdef.o igotdata.o inarea.o inlayer.o isitob.o readcntl.o \
            setarea.o setfcst.o setlevel.o setmodel.o setobtyp.o stclst.o stldsp.o \
            stlstr.o stnumb.o stread.o strmbl.o strxbl.o \
            w3fb06.o w3fb11.o w3fb71.o gridtobs.o
CMD=        verf_gridtobs_gridtobs_ens

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS)$(LDFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
