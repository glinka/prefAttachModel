toCompile = prefAttachModel.o main.o

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x -O3 #-fPIC

all: prefAttachModel

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

prefAttachModel: $(toCompile)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) *.o 
