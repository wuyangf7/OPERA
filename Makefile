
# -----------------------------------------------------------------
#   Makefile for OPERA 
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------

# Directory of the target
OUTPUT = opera_Linux

# Compiler
CXX = g++

# EIGEN library
EIGEN_PATH = $(EIGEN3_INCLUDE_DIR)
BOOST = /home/uqywu16/90days/program/lib/boost_1_61_0

# Intel MKL library
#MKL_PATH = /opt/intel/mkl

# Compiler flags
CXXFLAGS = -w -O3 -m64 -fopenmp -I $(BOOST) -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -std=c++11 
#LIB += -static -lz -Wl,-lm -ldl
LIB += -static -lz -lgomp

HDR += CommFunc.h \
	   cdflib.h \
	   dcdflib.h \
           SMR.h \
	   ipmpar.h \
           StatFunc.h \
           StrFunc.h \
            SMR_data.h \
            SMR_plot.h \
            SMR_data_p1.h \
            SMR_data_p2.h \
	   stat.hpp \
	   bfile.hpp \
           SMR_data_p3.h 
SRC = SMR.cpp \
           CommFunc.cpp \
           SMR_data.cpp \
	   dcdflib.cpp \
           StatFunc.cpp \
           StrFunc.cpp	\
           SMR_plot.cpp \
           SMR_data_p1.cpp \
           SMR_data_p2.cpp \
	   stat.cpp \
	   bfile.cpp \
           SMR_data_p3.cpp 
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f *.o
