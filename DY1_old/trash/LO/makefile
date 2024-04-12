FC = gfortran
FFLAGS   = -O
LDFLAGS = -L./../ -lfile
export LHFLAGS=$(shell lhapdf-config --ldflags)
OBJECTS =alphaS.o misc.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o
PROG = runLO
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
clean: 
	rm -f *.o core *~ *.x $(PROG)
