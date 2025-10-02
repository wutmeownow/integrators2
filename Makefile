# here we access the root configuration, include files, and libraries
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMathMore
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  = $(ROOTCFLAGS) -I$(ODELIB) -g #-Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS)
GXX	   = g++ $(CXXFLAGS)


a11: sobol 5630start 3630start

ndcrescent_plotting: ndcrescent_plotting.cpp
	$(GXX) $(CXXFLAGS) -o ndcrescent_plotting ndcrescent_plotting.cpp $(LDFLAGS)

ndcrescent: ndcrescent.cpp
	$(GXX) $(CXXFLAGS) -o ndcrescent ndcrescent.cpp $(LDFLAGS)

sobol: sobol.cpp 
	$(GXX) $(CXXFLAGS) -o sobol sobol.cpp $(LDFLAGS)

5630start: 5630start.cpp
	$(GXX) $(CXXFLAGS) -o 5630start 5630start.cpp $(LDFLAGS)

3630start: 3630start.cpp
	$(GXX) $(CXXFLAGS) -o 3630start 3630start.cpp $(LDFLAGS)

clean:
	rm -f sobol 56xxstart 36xxstart
	rm -f *~ *.d *.so *.pcm 

cleanall: clean
	rm -f *png *pdf 
