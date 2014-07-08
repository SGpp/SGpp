export SGPP=$(HOME)/projs/sg/SGpp/branches/sg-online-learning

SRCDIR=$(SGPP)/src/sgpp
OBJDIR=$(SGPP)/tmp
LIBDIR=$(SGPP)/lib/sgpp
INCLDIR=$(SGPP)/src/sgpp
GA_INCLDIR=$(HOME)/opt/include

INCLUDE_FLAGS=-I$(INCLDIR) -I$(GA_INCLDIR) $(GSL_INC)

TARGET=main

CXXFLAGS = -g -O0 -std=c++0x -Wno-deprecated
LDFLAGS  = -g -L/usr/lib
CXX=g++
LDLIBS=-lsgppbase -lsgppdatadriven -lsgpppde -lsgppsolver -lsgppcombigrid

default: $(TARGET)

$(TARGET): $(TARGET).cpp PerformanceTest.cpp test/TestData.cpp test/FriedmanTestData.cpp
	$(CXX) $(CXXFLAGS) -L./ $(LDFLAGS) -L$(LIBDIR) $(INCLUDE_FLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.so *.o $(TARGET) test/*.so test/*.o
