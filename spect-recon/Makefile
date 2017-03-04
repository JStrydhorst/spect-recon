MAKE = make
CXX = g++
CXXFLAGS = -c -std=c++11 -Wall -O3
LFLAGS = -Wall
INC_DIR = ../lib
LIB_DIR = ../lib
OBJS = spect_osem4d.o AttenMap.o Projection.o Reconstruction.o
LIBS = $(addprefix $(LIB_DIR)/, fft.o fft8g.o HDLListObj.o DataElement.o DicomObj.o RootDicomObj.o ProjTable.o Pinhole.o)

.PHONY : all libs clean

all : spect_recon

$(LIBS) :
	$(MAKE) -C $(LIB_DIR)

spect_osem4d.o : spect_osem4d.cpp Reconstruction.h Projection.h version.h
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $< 

Reconstruction.o : Projection.h AttenMap.h

spect_recon : $(OBJS) $(LIBS)
	$(CXX) $(LFLAGS) $(OBJS) $(LIBS) -o spect_recon

%.o : %.cpp %.h 
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $< 

clean :
	@rm -f *.o
	@$(MAKE) -C $(LIB_DIR) clean
	
