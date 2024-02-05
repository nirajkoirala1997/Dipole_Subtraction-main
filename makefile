FC = gfortran
FFLAGS   = -O
#LDFLAGS = -L./ -ldhelas3
LIBS = -lcuba
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
LIB = 
OBJECTS = cuba.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
PROG = run.x
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS)
clean: 
	rm -f *.o core *~ *.x
#$(PROG): *.inc
#$(PROG).o: $(PROG).f param_card.dat
