MODEL_SRCS=main.cc prefAttachModel.cc pamCPI.cc calcGraphProps.cc fitCurves.cc util_fns.cc util_fns.cc
NG_SRCS=coarse_newton_main.cc newton_wrapper.cc newton.cc gmres.cc prefAttachModel.cc pamCPI.cc calcGraphProps.cc fitCurves.cc util_fns.cc 
GE_SRCS=graph_embedding_main.cc util_fns.cc prefAttachModel.cc eigen_solvers.cc
MODEL_OBJECTS=$(MODEL_SRCS:.cc=.o)
NG_OBJECTS=$(NG_SRCS:.cc=.o)
GE_OBJECTS=$(GE_SRCS:.cc=.o)

# CXX = g++
CXXFLAGS = -I ~/build/Eigen -debug full -std=c++0x -mkl -gxx-name=/usr/bin/g++ -O3 -traceback #-openmp # #/home/oakridge/holiday/build/bin/g++

# CXX = g++
# CXXFLAGS = -g -Wall -Wno-sign-compare -std=c++0x #-O3

CXX = mpic++
# CXXFLAGS = -g -Wall -Wno-sign-compare -std=c++0x #-O3

all: graph_embedding # pref_attach # coarse_ng 

%.o: %.c
	$(CXX) -c $<  $(CXXFLAGS)

pref_attach: $(MODEL_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

coarse_ng: $(NG_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

graph_embedding: $(GE_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

depend: .depend

.depend: $(MODEL_SRCS) $(NG_SRCS) $(GE_SRCS)
	rm -f ./.depend
	$(CXX) -MM -MT $^ $(CXXFLAGS) > ./.depend

clean:
	$(RM) *.o 

include .depend
