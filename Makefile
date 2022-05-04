ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTINCLUDE  := -I$(shell root-config --incdir)

exe = MCEGENpiN_radcorr

all: $(exe)

MCEGENpiN_radcorr:
	$(CXX) -O3 $(ROOTINCLUDE) $(ROOTCFLAGS) -o $(exe) *.cpp $(ROOTLIBS)

clean:
	rm -rf $(exe)

