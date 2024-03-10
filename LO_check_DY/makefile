FC = gfortran
FFLAGS   = -O
LDFLAGS = -lc -L./../ -lfile
LIBS = -lcuba
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
OBJECTS =alphaS.o misc.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o
##OBJECTS =alphaS.o misc.o lum.o pdf_lha.o vsup.o integrand.o kin2.o
PROG = runLO
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
lum.o: ../lum.f
	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
vsup.o: ../vsup.f
	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
pdf_lha.o: ../pdf_lha.f
	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
#cuba.o: ../cuba.f
#	$(FC) $(FFLAGS) -c ../cuba.f -o cuba.o
alphaS.o: ../alphaS.f
	$(FC) $(FFLAGS) -c ../alphaS.f -o alphaS.o
phasespace.o: ../phasespace.f
	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
misc.o: ../misc.f
	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
clean: 
	rm -f *.o core *~ *.x $(PROG)
$(PROG).o: $(PROG).f param_card.dat
