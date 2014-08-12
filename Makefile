SRCS=main.cc prefAttachModel.cc pamCPI.cc calcGraphProps.cc fitCurves.cc util_fns.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = icpc 

CXXFLAGS = -I/home/holiday/build/Eigen -debug full -std=c++0x -O3 -mkl -gxx-name=/usr/bin/g++ #-openmp # #/home/oakridge/holiday/build/bin/g++
# CXX = g++

# CXXFLAGS = -g -Wall -Wno-sign-compare -std=c++0x #-O3

all: pref_attach

%.o: %.c
	$(CXX) -c $<  $(CXXFLAGS)

pref_attach: $(OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) -MM -MT $^ $(CXXFLAGS) > ./.depend

clean:
	$(RM) *.o 

include .depend
