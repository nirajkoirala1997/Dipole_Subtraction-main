CXX = g++
F77=gfortran
FFLAGS= -c
PREFIXG=/home/niraj/Vaibhav/installations

#LDFLAGS = -I$(PREFIXG)/include -L$(PREFIXG)/lib -lginac -L$(PREFIXG)/lib -lcln -lstdc++
#CXXFLAGS = -DDIGITS=30 -I$(PREFIXG)/include -L$(PREFIXG)/lib -lginac -std=c++11
LDFLAGS =  -lginac -lcln -lstdc++
CXXFLAGS = -DDIGITS=30  -lginac -std=c++11

#load module compilers/gcc/9.3.0

export LHFLAGS=$(shell lhapdf-config --prefix)/include
export LHAPDF=$(shell lhapdf-config --prefix)/lib
export handyGL=/home/niraj/Vaibhav/installations/lib
export handyGI=/home/niraj/Vaibhav/installations/include

obj=$(wildcard *.f)

run.x: evaluate_function.o cpp_wrapper.o  $(obj)
	$(F77) evaluate_function.o cpp_wrapper.o -o $@ $(obj) $(LDFLAGS) -L$(handyGL) -lhandyg  -I$(handyGI) -L$(LHAPDF) -lLHAPDF -I$(LHFLAGS)
%.o:%.cpp
	$(CXX) -c $(CXXFLAGS)   $< -o $@

clean: 
	rm -f *.o core *~ *.x
