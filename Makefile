THIS_DIR = ./lib
LOCAL_DIR = $(HOME)/local

XX = gcc
MARCH = native
STD = c++11 
# GCC flags, including many useful warnings
# (13.05.2018 @ 13:50) Options updated after reviewing the following GCC pages
# "Intel 386 and AMD x86-64 Options"
# "Options to Request or Suppress Warnings"
# "Options That Control Optimization"

# -pedantic-errors is like -Wpedantic, except the warnings become errors
STABILITY_FLAGS = -pedantic-errors -fno-common -mfpmath=sse -mieee-fp #sse flag to avoid weird x87 registers (see https://gcc.gnu.org/wiki/FloatingPointMath)
STABILITY_WARNINGS = -Wall -Wextra -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wformat -Wmissing-declarations -Wredundant-decls -Wmissing-field-initializers -Wlogical-op -Wundef -Wfloat-equal -Wdouble-promotion -Wstrict-aliasing # -Wunsafe-loop-optimizations -Wuseless-cast
PERFORMANCE_FLAGS = -O2 -march=$(MARCH) -Winline -Wdisabled-optimization -Wpadded -ftree-vectorize -funsafe-loop-optimizations -Wvector-operation-performance # vectorize is the only thing from O3 that we want
# Append processor dependent 
PERFORMANCE_FLAGS += $(shell sh getSSE_AVX.sh)
BUILD_LIB_FLAGS = -fPIC
# 
CXXFLAGS = -std=$(STD) $(STABILITY_WARNINGS) $(PERFORMANCE_FLAGS) $(BUILD_LIB_FLAGS) -g #-fopt-info-vec-optimized

# The directory structure of pqRand
INCLUDE = ./include
SOURCE = ./source
EXAMPLES = ./examples

INC_FLAGS = -I $(INCLUDE) -I $(LOCAL_DIR)/include
LOCAL_LIBS = -lpqr
EXTERN_LIB_FLAGS = -lstdc++ -lm -lQtCore -L $(LOCAL_DIR)/lib $(LOCAL_LIBS)
LIB_FLAGS = $(EXTERN_LIB_FLAGS) -L $(THIS_DIR) -lkdp

EXAMPLES_CPP = $(wildcard testing/*.cpp)
EXAMPLES_X = $(patsubst %.cpp, %.x, $(EXAMPLES_CPP))

FILENAMES = kdpVectors kdpHistogram
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
