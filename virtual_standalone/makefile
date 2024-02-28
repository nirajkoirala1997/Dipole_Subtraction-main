FC = gfortran
FFLAGS   = -O
LDFLAGS = -L./ -ldhelas3 -lc -L./../ -lfile
LIBS = -lcuba
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
OBJECTS = misc.o Iterm.o matrixLO1.o couplings.o cmatrix.o cuba.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o
PROG = runVir
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
clean: 
	rm -f *.o core *~ *.x $(PROG)
$(PROG): *.inc
$(PROG).o: $(PROG).f param_card.dat
