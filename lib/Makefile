CC = gcc
CXX = g++
CFLAGS = -c -Wall
CXXFLAGS = -c -Wall -std=c++11

OBJS := fft.o fft8g.o HDLListObj.o DataElement.o DicomObj.o RootDicomObj.o ProjTable.o Pinhole.o

.PHONY : all clean

all : $(OBJS)

DataElement.o : HDLListObj.h DicomObj.h RootDicomObj.h
DicomObj.o : HDLListObj.h DataElement.h
RootDicomObj.o : HDLListObj.h DicomObj.h DataElement.h

%.o : %.c
	$(CC) $(CFLAGS) $<

%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) $<

clean :
	@rm -f *.o

