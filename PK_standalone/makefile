FC = gfortran#g77
#FC = f77#gfortran#g77
FFLAGS   = -O 
LDFLAGS = -lc -L./../ -lfile
MYFLAG = #-L./home/niraj/1TB-Disc/Workspace-IITG/Packages/Source/mylibrary -lfile
#export LHFLAGS = $(shell lhapdf-config --ldflags)
LHFLAGS = $(shell lhapdf-config --ldflags)
PROG    = runPK
OBJECTS4  = alphaS.o PK_new.o pdf_lha.o lum.o integrand.o phasespace.o vsup.o
	
$(PROG): $(OBJECTS4) main.f
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS4) $(LDFLAGS) $(LHFLAGS) $(MYFLAG) 
#	rm -f *.o
clean:
	rm -f $(OBJECTS4)
	rm -f $(PROG).o $(PROG) 

distclean:
	make clean
	rm -f  $(PROG)

$(PROG).o :  $(PROG).f param_card.dat
cmatrix.o: *.inc

