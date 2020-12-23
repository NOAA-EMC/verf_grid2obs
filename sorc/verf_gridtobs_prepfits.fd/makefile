SHELL =      /bin/sh

FC =         ifort

LIBS    = $(G2_LIB4) $(W3NCO_LIB4) $(W3EMC_LIB4) $(BUFR_LIB4_DA) $(BACIO_LIB4) $(JASPER_LIB) $(PNG_LIB) $(Z_LIB)


INC     = $(G2_INC4)

FFLAGS =      -O2 -I$(INC) -mcmodel=medium -shared-intel -g -traceback -debug extended

#LDFLAGS =    -shared-intel 

SRCS=        anl.f grid3d.f surfce.f gridef.f observ.f backgv.f guesfc.f guess.f itotal.f obstrs.f weights.f debug.f counts.f vdates.f vrtfac.f guser.f adddate.f calcape.f calpbl.f calpw.f  getanl.f getbak.f getfct.f getprof_upr.f htindex.f ireadpb.f \
             obscape.f obspbl.f  calshear.f raddate.f setterp.f ufbxy3.f  windchill.f  haines.f   prepfits.f
OBJS=        anl.o grid3d.o surfce.o gridef.o observ.o backgv.o guesfc.o guess.o itotal.o obstrs.o  weights.o debug.o counts.o vdates.o vrtfac.o guser.o adddate.o calcape.o calpbl.o calpw.o  getanl.o getbak.o getfct.o getprof_upr.o htindex.o ireadpb.o \
             obscape.o obspbl.o  calshear.o raddate.o setterp.o ufbxy3.o windchill.o haines.o prepfits.o
CMD=            verf_gridtobs_prepfits

all:	$(CMD)

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS)$(LDFLAGS) -g  -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -rf *.o *.mod *.lst
