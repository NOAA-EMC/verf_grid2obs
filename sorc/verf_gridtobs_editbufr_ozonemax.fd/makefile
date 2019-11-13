LIBS	=  $(W3EMC_LIB4) $(W3NCO_LIB4) $(BUFR_LIB4) $(BACIO_LIB4)
OBJS	= gtgdef.o ireadpb.o raddate.o w3fi63.o w3fi71.o editbufr.o
FC	= ifort
FFLAGS	= -heap-arrays
CMD	= verf_gridtobs_editbufr_ozonemax

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	rm -f *.o $(CMD)





