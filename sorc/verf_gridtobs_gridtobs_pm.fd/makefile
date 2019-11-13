OBJS	= gtgdef.o igotdata.o inarea.o inlayer.o isitob.o raddate.o readcntl.o \
	setarea.o setfcst.o setlevel.o setmodel.o setobtyp.o stclst.o stldsp.o \
	stlstr.o stnumb.o stread.o strmbl.o strxbl.o ufbin3.o w3fb06.o w3fb11.o \
	w3fi63.o w3fi71.o gets1loc.o gridtobs.o

SHELL=          /bin/sh
LIBS=            $(W3EMC_LIB4) $(W3NCO_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
FC=             ifort
FFLAGS=         -heap-arrays -i4
CMD	= verf_gridtobs_gridtobs_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	rm -f *.o $(CMD)

