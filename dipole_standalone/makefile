FC = gfortran
FFLAGS   = -O
LDFLAGS = -L./../ -lfile
LIBS = -lcuba
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
LIB = 
#OBJECTS = Iterm.o matrixLO1.o couplings.o cmatrix.o PK.o cuba.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
OBJECTS = alphaS.o cuba.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
PROG = runDipole
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
clean: 
	rm -f *.o core *~ *.x
$(PROG): 
$(PROG).o: $(PROG).f param_card.dat
