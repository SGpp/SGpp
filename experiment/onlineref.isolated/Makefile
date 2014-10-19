ifndef SGPP
  $(error "SGPP not set.")
endif

SRCDIR=$(SGPP)/src/sgpp
OBJDIR=$(SGPP)/tmp
LIBDIR=$(SGPP)/lib/sgpp
INCLDIR=$(SGPP)/src/sgpp
GA_INCLDIR=$(HOME)/opt/include

INCLUDE_FLAGS=-I$(INCLDIR) -I$(GA_INCLDIR) $(GSL_INC)

TARGET=main

CXXFLAGS = -g -O0 -std=c++0x -Wno-deprecated
#CXXFLAGS = -O3 -std=c++0x -Wno-deprecated
LDFLAGS  = -g -L/usr/lib
LDFLAGS  = -L/usr/lib
CXX=g++
LDLIBS=-lsgppbase -lsgppdatadriven -lsgpppde -lsgppsolver -lsgppcombigrid

default: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -L./ $(LDFLAGS) -L$(LIBDIR) $(INCLUDE_FLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.so *.o $(TARGET) test/*.so test/*.o
