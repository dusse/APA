
HDF5_PATH=.../libs/hdf5/
MPI_PATH=.../libs/openmpi/
PYTHON27_INC=.../include/python2.7/
PYTHON27_LIB=.../lib/python2.7/


LIBS=-lpython2.7 -lhdf5
INCLUDES=-I$(PYTHON27_INC) -I$(HDF5_PATH)/include/ -I$(MPI_PATH)/include/

DSRC = ./src
DEXE = ./

LD_LIBRARY_PATH=$(HDF5_PATH)/lib/:$(PYTHON27_LIB)
export LIBRARY_PATH=$LIBRARY_PATH:$(LD_LIBRARY_PATH)


CXX = $(MPI_PATH)/bin/mpicxx
CXXFLAGS  = -DLOG -O3 -Wall -c -std=c++11 -Wno-sign-compare -Wno-unused-variable

_SRCS =  $(DSRC)/core/SimulationManager.cpp \
               $(DSRC)/grid/GridManager.cpp \
               $(DSRC)/input/Loader.cpp \
               $(DSRC)/misc/Logger.cpp \
               $(DSRC)/misc/Misc.cpp \
               $(DSRC)/output/Writer.cpp \
               $(DSRC)/common/variables/VectorVar.cpp \
               $(DSRC)/solvers/EquationSolver.cpp \
               $(DSRC)/solvers/ModelInitializer.cpp \
               $(DSRC)/APA.cpp \

_OBJS            = $(_SRCS:.cpp=.o)

_EXEN            = $(DEXE)/apa.exe

all : $(_EXEN)


$(_EXEN) : $(_OBJS)
	@echo 'Building target: $@'
	$(CXX) -o $@ $^  $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o : %.cpp
	$(CXX) $(INCLUDES) -o $@ $< $(CXXFLAGS) 


clean :
	rm -f $(_OBJS)


