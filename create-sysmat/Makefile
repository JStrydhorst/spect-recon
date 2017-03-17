MAKE = make
CXX = g++
CXXFLAGS = -c -std=c++11 -Wall -msse4.2 -O3
LFLAGS = -Wall
INC_DIR = ../lib
LIB_DIR = ../lib
OBJS = create-sysmat.o 
LIBS = $(addprefix $(LIB_DIR)/, ProjTable.o Pinhole.o)

.PHONY : all libs clean

all : create-sysmat 

$(LIBS) :
	$(MAKE) -C $(LIB_DIR)

create-sysmat : $(OBJS) $(LIBS)
	$(CXX) $(LFLAGS) $(OBJS) $(LIBS) -o create-sysmat

%.o : %.cpp %.h 
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $< 

clean :
	@rm -f *.o
	@$(MAKE) -C $(LIB_DIR) clean
	
