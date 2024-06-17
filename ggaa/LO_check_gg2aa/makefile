# Compiler and flags
FC = gfortran
FFLAGS = -O
LDFLAGS = -lc -L./../ -lfile
export LHFLAGS = $(shell lhapdf-config --ldflags)

# Program and object files
PROG = runLO
OBJECTS = mat.amp.o cuts.lo2.o cone.o printframe.o misc.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o couplings.o

# Main target
$(PROG): $(OBJECTS) main.f
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS)

# Individual object file rules
lum.o: ../lum.f
	$(FC) $(FFLAGS) -c ../lum.f -o lum.o

vsup.o: ../vsup.f
	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o

pdf_lha.o: ../pdf_lha.f
	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o

phasespace.o: ../phasespace.f
	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o

misc.o: ../misc.f
	$(FC) $(FFLAGS) -c ../misc.f -o misc.o

printframe.o: ../printframe.f
	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o

mat.amp.o: ../mat.amp.f
	$(FC) $(FFLAGS) -c ../mat.amp.f -o mat.amp.o

couplings.o: ../slicing_files/couplings.f
	$(FC) $(FFLAGS) -c ../slicing_files/couplings.f -o couplings.o

cuts.lo2.o: ../slicing_files/cuts.lo2.f
	$(FC) $(FFLAGS) -c ../slicing_files/cuts.lo2.f -o cuts.lo2.o

cone.o: ../slicing_files/cone.f
	$(FC) $(FFLAGS) -c ../slicing_files/cone.f -o cone.o

# Clean target
clean:
	rm -f *.o core *~ *.x $(PROG)

