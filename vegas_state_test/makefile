# Compiler and flags
FC = gfortran
FFLAGS = -O -w

# Directories
OBJDIR = objects
SRCDIR = .
BINDIR = .

# Program and object files
PROG = $(BINDIR)/runLO
OBJECTS = $(OBJDIR)/vsup.o

# Ensure the object directory exists
$(OBJDIR)/%.o: $(SRCDIR)/%.f
	mkdir -p $(OBJDIR)
	$(FC) -c $(FFLAGS) $< -o $@

# Main target: Build executable from object files
$(PROG): $(OBJECTS) $(SRCDIR)/main.f
	$(FC) $(FFLAGS) $(OBJECTS) $(SRCDIR)/main.f -o $(PROG)

# Clean target
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)

