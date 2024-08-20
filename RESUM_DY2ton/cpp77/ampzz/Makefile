CXX=g++
# if DIGITS is undefined, double precision is used
# set -DDIGITS=30 to use 30 digits for internal evaluations
# this might be useful in extreme kinematic regimes
PREFIX=<PATH-TO-GiNaC-installation>

CXXFLAGS=-DDIGITS=30 -O1 -I$(PREFIX)/include 
LDFLAGS=-L$(PREFIX)/lib -lginac -lcln

OBJ=ampvv.o \
 ampvv01.o ampvv02.o ampvv03.o ampvv04.o ampvv05.o \
 ampvv06.o ampvv07.o ampvv08.o ampvv09.o ampvv10.o \
 ampvv11.o ampvv12.o ampvv13.o ampvv14.o ampvv15.o \
 ampvv16.o ampvv17.o ampvv18.o ampvv19.o ampvv20.o \
 ampvv21.o ampvv22.o ampvv23.o ampvv24.o ampvv25.o \
 ampvv26.o ampvv27.o ampvv28.o ampvv29.o ampvv30.o \
 ampvv31.o ampvv32.o ampvv33.o ampvv34.o ampvv35.o \
 ampvv36.o ampvv37.o ampvv38.o ampvv39.o ampvv40.o \
 ampvv41.o ampvv42.o ampvv43.o ampvv44.o ampvv45.o \
 ampvv46.o ampvv47.o ampvv48.o ampvv49.o ampvv50.o \
 ampvv51.o ampvv52.o ampvv53.o ampvv54.o ampvv55.o \
 ampvv56.o ampvv57.o ampvv58.o ampvv59.o ampvv60.o \
 ampvv61.o ampvv62.o ampvv63.o ampvv64.o ampvv65.o \
 ampvv66.o ampvv67.o ampvv68.o ampvv69.o ampvv70.o \
 ampvv71.o ampvv72.o ampvv73.o ampvv74.o ampvv75.o \
 ampvv76.o ampvv77.o ampvv78.o ampvv79.o ampvv80.o \
 ampvv81.o ampvv82.o ampvv83.o ampvv84.o ampvv85.o \
 ampvv86.o ampvv87.o ampvv88.o ampvv89.o ampvv90.o \
 ampvv91.o ampvv92.o ampvv93.o ME2_ZZ.o

all: test
	ar -r libZZ.a *.o

test: test.o $(OBJ)
	$(CXX) test.o $(OBJ) -o $@ $(LDFLAGS)

libZZ: ME2_ZZ.o $(OBJ)
#	$(CXX) test.o $(OBJ) -o $@ $(LDFLAGS)

ampvv.cpp: ampvva.cpp ampvvb.cpp ampvvc.cpp
	cat ampvva.cpp ampvvb.cpp ampvvc.cpp > ampvv.cpp

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -f ampvv.cpp *.o *~ test

.PHONY: clean
