SRCS=main.cc prefAttachModel.cc pamCPI.cc calcGraphProps.cc fitCurves.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -Wno-sign-compare -std=c++0x -O3

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
