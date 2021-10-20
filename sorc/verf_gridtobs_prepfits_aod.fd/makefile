SHELL=          /bin/sh
LIBS =   $(G2_LIB4) $(W3EMC_LIB4) $(W3NCO_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)\
        -L/usrx/local/64bit.cirrus/lib $(JASPER_LIB) $(PNG_LIB) $(Z_LIB)
FC=             ftn
INC =           $(G2_INC4)
FFLAGS=         -heap-arrays -i4 -O0 -I$(INC)
LDFLAGS=       
SRCS=        grid3d.f surfce.f gridef.f observ.f backgv.f guesfc.f guess.f itotal.f obstrs.f weights.f debug.f counts.f vdates.f vrtfac.f guser.f satimg.f adddate.f calcape.f calpbl.f calpw.f  getbak.f getfct.f getprof_upr.f htindex.f ireadpb.f \
             obscape.f obspbl.f  raddate.f setterp.f ufbxy3.f ufbxy3_goes.f w3fa09a.f w3fb11a.f w3fb14.f w3fi71.f windchill.f  prepfits.f
OBJS=        grid3d.o surfce.o gridef.o observ.o backgv.o guesfc.o guess.o itotal.o obstrs.o  weights.o debug.o counts.o vdates.o vrtfac.o guser.o satimg.o adddate.o calcape.o calpbl.o calpw.o  getbak.o getfct.o getprof_upr.o htindex.o ireadpb.o \
             obscape.o obspbl.o  raddate.o setterp.o ufbxy3.o ufbxy3_goes.o w3fa09a.o w3fb11a.o w3fb14.o w3fi71.o windchill.o prepfits.o
CMD=            ../../exec/verf_gridtobs_prepfits_aod

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS)$(LDFLAGS) -g  -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.lst
