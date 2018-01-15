THIS_DIR = ./lib
LOCAL_DIR = $(HOME)/local

CXX = gcc
MARCH = native
STD = c++11 
# GCC flags, including many useful warnings
STABILITY_FLAGS = -pedantic-errors -fno-common -mfpmath=sse -mieee-fp #sse flag to avoid weird x87 registers (see https://gcc.gnu.org/wiki/FloatingPointMath)
STABILITY_WARNINGS = -Wall -Wextra -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wmissing-declarations -Wredundant-decls -Wmissing-field-initializers -Wlogical-op -Wunsafe-loop-optimizations -Wwrite-strings -Wundef -Wfloat-equal
PERFORMANCE_FLAGS = -O2 -march=$(MARCH) -msse4 -mavx2 -Winline -Wdisabled-optimization -Wpadded -ftree-vectorize # vectorize is the only thing from O3 that we want
BUILD_LIB_FLAGS = -fPIC
# 
CXXFLAGS = -std=$(STD) $(STABILITY_WARNINGS) $(PERFORMANCE_FLAGS) $(BUILD_LIB_FLAGS) -g #-fopt-info-vec-optimized

# The directory structure of pqRand
INCLUDE = ./include
SOURCE = ./source
EXAMPLES = ./examples

INC_FLAGS = -I $(INCLUDE) -I $(LOCAL_DIR)/include
LOCAL_LIBS =
EXTERN_LIB_FLAGS = -lstdc++ -lm -L $(LOCAL_DIR)/lib $(LOCAL_LIBS)
LIB_FLAGS = $(EXTERN_LIB_FLAGS) -L $(THIS_DIR) -lkdp

EXAMPLES_CPP = $(wildcard examples/*.cpp)
EXAMPLES_X = $(patsubst %.cpp, %.x, $(EXAMPLES_CPP))

FILENAMES = kdpVectors
OBJS = $(addsuffix .o, $(addprefix $(SOURCE)/, $(FILENAMES)))

all : lib $(EXAMPLES_X)

lib : $(THIS_DIR)/libkdp.so

$(THIS_DIR)/libkdp.so: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(EXTERN_LIB_FLAGSs) -o $@

%.x : %.cpp $(THIS_DIR)/libkdp.so
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) $(LIB_FLAGS) $*.cpp -o $@
	
%.o : %.cpp 
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) $(EXTERN_LIB_FLAGS) $*.cpp -c -o $*.o
	
.PHONY: clean

clean:
	rm -f $(SOURCE)/*.o
	rm -f $(THIS_DIR)/libkdp.so
	rm -f examples/*.x
	rm -f testing/*.x
