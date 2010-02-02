
F77    = ifort 
FFLAGS = 
LINK   = ifort
LFLAGS = -static-libcxa

OBJS   = fit3pbnd.o dlsei.o dlsi.o dwnnls.o dh12.o dhfti.o dlpdp.o pcmach.o \
	 math.o xerror.o

fit3pbnd.exe : $(OBJS)
	$(LINK) $(LFLAGS) -o fit3pbnd.exe $(OBJS)
	cp fit3pbnd.exe ../run

.f.o:
	$(F77) -c $(FFLAGS) $<

clean : 
	rm -f *.exe *.a *.o
