# Compiler and Flags
CXX = g++
F77 = gfortran
FFLAGS = -c
PREFIXG=/home/niraj/Vaibhav/installations
# Paths

LDFLAGS =  -lginac -lcln -lstdc++
CXXFLAGS = -DDIGITS=30  -lginac -std=c++11

export LHFLAGS=$(shell lhapdf-config --prefix)/include
export LHAPDF=$(shell lhapdf-config --prefix)/lib
export handyGL=/home/niraj/Vaibhav/installations/lib
export handyGI=/home/niraj/Vaibhav/installations/include

# Object Files
FSRCS = $(wildcard *.f)
FOBJS = $(FSRCS:.f=.o)

obj=$(wildcard *.f)

run.x: evaluate_function.o cpp_wrapper.o  $(obj)
	$(F77) evaluate_function.o cpp_wrapper.o -o $@ $(obj) $(LDFLAGS) -L$(handyGL) -lhandyg  -I$(handyGI) -L$(LHAPDF) -lLHAPDF -I$(LHFLAGS)
%.o:%.cpp
	$(CXX) -c $(CXXFLAGS)   $< -o $@

clean: 
	rm -f *.o core *~ *.x


# Targets
#run.x: evaluate_function.o cpp_wrapper.o $(FOBJS)
#	$(F77) -o $@ evaluate_function.o cpp_wrapper.o $(FOBJS) $(LDFLAGS) -L$(handyGL) -lhandyg  -I$(handyGI) -L$(LHAPDF) -lLHAPDF -I$(LHFLAGS)
#
#%.o: %.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@
#
#%.o: %.f
#	$(F77) $(FFLAGS) $< -o $@
#
#clean:
#	rm -f *.o run.x
#


#
#CXX = g++
#F77=gfortran
#FFLAGS= -c
#Pathinstall = /home/chinmoy/Work/Packages/Install
#
#CXXFLAGS = -DDIGITS=30 -lginac -lcln
#LDFLAGS =  -lginac -lcln -lstdc++ 
##LDFLAGS = -I$(Pathinstall)/include  -L$(Pathinstall)/lib -lginac -lcln -lstdc++
#LDFLAGS = 
#
#export LHFLAGS=$(shell lhapdf-config --prefix)/include
#export LHAPDF=$(shell lhapdf-config --prefix)/lib
##export handyGL=$(Pathinstall)/lib
##export handyGI=$(Pathinstall)/include
##export handyGL=/home/vaibhav/work/Packages/handyg-master/handyGnew/handyg/install/lib
##export handyGI=/home/vaibhav/work/Packages/handyg-master/handyGnew/handyg/install/include
#
#obj=$(wildcard *.f)
#
#run.x: evaluate_function.o cpp_wrapper.o  $(obj)
#	$(F77) evaluate_function.o cpp_wrapper.o -o $@ $(obj) $(LDFLAGS) -L$(LHAPDF) -lLHAPDF -I$(LHFLAGS)$  
##	$(F77) evaluate_function.o cpp_wrapper.o -o $@ $(obj) $(LDFLAGS) -L$(handyGL) -lhandyg  -I$(handyGI) -L$(LHAPDF) -lLHAPDF -I$(LHFLAGS)$  
#%.o:%.cpp
#	$(CXX) -c $(CXXFLAGS) -std=c++14 $< -o $@
#
#clean: 
#	rm -f *.o core *~ *.x
