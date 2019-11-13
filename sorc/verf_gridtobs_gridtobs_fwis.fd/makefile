SHELL=          /bin/sh
LIBS =        $(G2_LIB4) $(JASPER_LIB) $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIB4) $(BACIO_LIB4) $(PNG_LIB) $(Z_LIB)
FC=             ifort
INC =        $(G2_INC4)
FFLAGS =     -i4 -g -traceback -I$(INC)
LDFLAGS =    -shared-intel
SRCS=       gtgdef.f igotdata.f inarea.f inlayer.f isitob.f readcntl.f \
            setarea.f setfcst.f setlevel.f setmodel.f setobtyp.f stclst.f stldsp.f \
            stlstr.f stnumb.f stread.f streadfho.f streal.f strmbl.f strxbl.f \
            maskdef.f gridtobs.f
OBJS=       gtgdef.o igotdata.o inarea.o inlayer.o isitob.o readcntl.o \
            setarea.o setfcst.o setlevel.o setmodel.o setobtyp.o stclst.o stldsp.o \
            stlstr.o stnumb.o stread.o streadfho.o streal.o strmbl.o strxbl.o \
            maskdef.o gridtobs.o
CMD=        verf_gridtobs_gridtobs_fwis

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS)$(LDFLAGS) -g  -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
