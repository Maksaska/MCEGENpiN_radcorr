CC=g++
CFLAGS=-c -Wall
LDFLAGS=-I$(shell root-config --incdir)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
SOURCES=App.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MCEGENpiN_radcorr

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(ROOTCFLAGS) $(OBJECTS) -o $@ $(ROOTLIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(LDFLAGS) $(ROOTCFLAGS) $< -o $@ $(ROOTLIBS)
	
clean:
	rm -rf *.o *.dat $(EXECUTABLE)
	
delete_data:
	rm *.dat *.jpeg
