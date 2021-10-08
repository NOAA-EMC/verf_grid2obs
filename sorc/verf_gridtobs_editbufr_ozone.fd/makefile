#!/bin/sh

OBJS    = gtgdef.o raddate.o w3fi63.o w3fi71.o editbufr.o
LIBS	=  $(W3EMC_LIB4) $(W3NCO_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
FC 	= ftn
FFLAGS 	= -heap-arrays 
CMD     = ../../exec/verf_gridtobs_editbufr_ozone

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS) 
clean:
	rm -f *.o $(CMD)
