FC = gfortran
FFLAGS   = -O
LDFLAGS = -L./../ -lfile
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
LIB = 
#OBJECTS = Iterm.o matrixLO1.o couplings.o cmatrix.o PK.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
OBJECTS = dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
PROG = runDipole
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
lum.o: ../lum.f
	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
vsup.o: ../vsup.f
	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
pdf_lha.o: ../pdf_lha.f
	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
misc.o: ../misc.f
	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
phasespace.o: ../phasespace.f
	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
clean: 
	rm -f *.o core *~ *.x $(PROG)
$(PROG): 
$(PROG).o: $(PROG).f param_card.dat
