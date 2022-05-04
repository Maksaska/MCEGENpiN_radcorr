ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTINCLUDE  := -I$(shell root-config --incdir)

exe = pi0p-pin-generator

all: $(exe)

pi0p-pin-generator:
	$(CXX) -O3 $(ROOTINCLUDE) $(ROOTCFLAGS) -o $(exe) *.cpp $(ROOTLIBS)

clean:
	rm -rf $(exe)

