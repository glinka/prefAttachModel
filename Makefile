MODEL_SRCS=main.cc prefAttachModel.cc pamCPI.cc calcGraphProps.cc fitCurves.cc custom_util_fns.cc
NG_SRCS=coarse_newton_main.cc newton_wrapper.cc newton.cc gmres.cc prefAttachModel.cc pamCPI.cc calcGraphProps.cc fitCurves.cc custom_util_fns.cc 
GE_SRCS=graph_embedding_main.cc custom_util_fns.cc prefAttachModel.cc eigen_solvers.cc calcGraphProps.cc
GEMOTIFS_SRCS=graph-embedding-motifs.cc custom_util_fns.cc prefAttachModel.cc eigen_solvers.cc calcGraphProps.cc
RHO_KAPPA_SRCS=kappa_rho_embedding_main.cc custom_util_fns.cc prefAttachModel.cc calcGraphProps.cc
MODEL_OBJECTS=$(MODEL_SRCS:.cc=.o)
NG_OBJECTS=$(NG_SRCS:.cc=.o)
GE_OBJECTS=$(GE_SRCS:.cc=.o)
GEMOTIFS_OBJECTS=$(GEMOTIFS_SRCS:.cc=.o)
RHO_KAPPA_OBJECTS=$(RHO_KAPPA_SRCS:.cc=.o)

# CXX = g++
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# may need to 'module purge intel/openmpi'
# then 'module load intel/openmpi'
# to place intel's 'mpic++' in its
# proper location in PATH
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CXXFLAGS = -I~/build/Eigen -I/home/oakridge/holiday/workspace/dmaps -I/home/oakridge/holiday/workspace/newton_gmres -I/home/oakridge/holiday/workspace/util_fns -I./igraph/include/igraph -L./igraph/lib -debug full -std=c++0x -mkl -gxx-name=/usr/bin/g++ -traceback -lutil_fns -leigensolvers -ligraph -openmp -O3 #/home/oakridge/holiday/build/bin/g++


# CXX = g++
# CXXFLAGS = -g -Wall -Wno-sign-compare -std=c++0x #-O3

CXX = mpic++
# CXXFLAGS = -g -Wall -Wno-sign-compare -std=c++0x #-O3

all: graph-embedding-motifs # graph_embedding # pref_attach coarse_ng rho_kappa_embedding

%.o: %.c
	$(CXX) -c $<  $(CXXFLAGS)

pref_attach: $(MODEL_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

coarse_ng: $(NG_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

graph_embedding: $(GE_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

rho_kappa_embedding: $(RHO_KAPPA_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

graph-embedding-motifs: $(GEMOTIFS_OBJECTS)
	$(CXX) -Wl,-rpath ./igraph/lib -I./igraph/include/igraph -L./igraph/lib -o $@ $^  $(CXXFLAGS) -ligraph

depend: .depend

.depend: $(MODEL_SRCS) $(NG_SRCS) $(GE_SRCS)
	rm -f ./.depend
	$(CXX) -MM -MT $^ $(CXXFLAGS) > ./.depend

clean:
	$(RM) *.o 

include .depend
