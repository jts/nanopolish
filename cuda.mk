#Make file options for CUDA support

NVCC ?= nvcc
CUDA_ROOT = /usr/local/cuda
CUDA_LIB ?= $(CUDA_ROOT)/lib64
CUDA_INCLUDE ?= $(CUDA_ROOT)/include
CURTFLAGS = -L$(CUDA_LIB) -lcudart
NVCCFLAGS ?= -std=c++11 -I. -I$(CUDA_INCLUDE) -O3 -use_fast_math --default-stream per-thread -restrict

CPPFLAGS += -I$(CUDA_INCLUDE)
CPPFLAGS += -DHAVE_CUDA=1

# Sub directories containing CUDA source code
SUBDIRS += src/cuda_kernels
# Find the source files by searching subdirectories
CU_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cu))
# Automatically generated object names
CU_OBJ = $(CU_SRC:.cu=.o)
CPP_OBJ += $(CU_OBJ)
LDFLAGS += $(CURTFLAGS)

.SUFFIXES: .cu

# Compile objects
.cu.o:
	$(NVCC) -o $@ -c $(NVCCFLAGS) $(CPPFLAGS) $<

