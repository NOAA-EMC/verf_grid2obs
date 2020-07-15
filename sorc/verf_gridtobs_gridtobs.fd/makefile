SHELL=          /bin/sh
LIBS=            $(W3EMC_LIB4) $(W3NCO_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
###LIBS=           -L/u/Jack.Woollen/lib -lw3_4 -lbufr_sc $(BACIO_LIB4)
FC=             ifort
FFLAGS=         -heap-arrays -i4 
LDFLAGS=       
SRCS=       gtgdef.f igotdata.f inarea.f inlayer.f isitob.f readcntl.f \
            setarea.f setfcst.f setlevel.f setmodel.f setobtyp.f stclst.f stldsp.f \
            stlstr.f stnumb.f stread.f streadfho.f streal.f strmbl.f strxbl.f \
            w3fa09a.f gridtobs.f
OBJS=       gtgdef.o igotdata.o inarea.o inlayer.o isitob.o readcntl.o \
            setarea.o setfcst.o setlevel.o setmodel.o setobtyp.o stclst.o stldsp.o \
            stlstr.o stnumb.o stread.o streadfho.o streal.o strmbl.o strxbl.o \
            w3fa09a.o gridtobs.o
CMD=        verf_gridtobs_gridtobs

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS)$(LDFLAGS) -g  -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
