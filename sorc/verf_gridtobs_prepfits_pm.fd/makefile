OBJS	= calcape.o adddate.o getbak.o getfct.o getprof.o ireadpb.o raddate.o \
	setterp.o ufbxy3.o w3fb13.o prepfits.o
SRCS	= calcape.f adddate.f getbak.f getfct.f getprof.f ireadpb.f raddate.f\
	setterp.f ufbxy3.f w3fb13.f prepfits.f
SHELL=          /bin/sh
LIBS =        $(G2_LIB4) $(JASPER_LIB) $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIBs) $(BACIO_LIB4) $(PNG_LIB) $(Z_LIB)
INC =        $(G2_INC4)
FFLAGS =     -i4 -g -traceback -I$(INC)
LDFLAGS =    -shared-intel
FC=             ifort
CMD	= verf_gridtobs_prepfits_pm

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	rm *.o
