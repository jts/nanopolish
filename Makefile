#

# Sub directories containing source code, except for the main programs
SUBDIRS := src src/hmm src/thirdparty src/thirdparty/scrappie src/common src/alignment src/pore_model src/io src/basemods

#
# Set libraries, paths, flags and options
#

#Basic flags every build needs
LIBS = -lz
CXXFLAGS ?= -g -O3
CXXFLAGS += -std=c++11 -fopenmp -fsigned-char -D_FILE_OFFSET_BITS=64 
CFLAGS ?= -O3 -std=c99 -fsigned-char -D_FILE_OFFSET_BITS=64
LDFLAGS ?=
CXX ?= g++
CC ?= gcc

ifeq ($(zstd),1)
LDFLAGS += -lzstd
endif

# Change the value of HDF5, EIGEN, or HTS below to any value to disable compilation of bundled code
HDF5 ?= install
EIGEN ?= install
HTS ?= install
MINIMAP2 ?= install
SLOW5LIB ?= install

# Check whether we're building on an ARM mac
ifneq (,$(findstring arm64-apple,$(MACHTYPE)))
	ARM_MAC=0
else
    ARM_MAC=1
    ARM=1
    SDK := $(shell xcrun --sdk macosx --show-sdk-path)
    MAC_FLAGS = -isysroot $(SDK) -L/usr/lib
endif

CXXFLAGS += $(MAC_FLAGS)
CFLAGS += $(MAC_FLAGS)

ifeq ($(ARM_MAC), 1)
    HDF5_VERSION ?= 1.13.0
else
    HDF5_VERSION ?= 1.8.14
endif

EIGEN_VERSION ?= 3.3.7

# Check operating system, OSX doesn't have -lrt
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBS += -lrt
endif

# Default to automatically installing hdf5
ifeq ($(HDF5), install)
    H5_LIB = ./lib/libhdf5.a
    H5_INCLUDE = -I./include
    LIBS += -ldl
else
    # Use system-wide hdf5
    H5_LIB =
    H5_INCLUDE ?=
    LIBS += -lhdf5
endif

# Default to automatically installing EIGEN
ifeq ($(EIGEN), install)
    EIGEN_CHECK = eigen/INSTALL
    EIGEN_INCLUDE = -I./eigen/
else
    # Use system-wide eigen
    EIGEN_CHECK =
    EIGEN_INCLUDE ?=
endif

# Default to build and link the libhts submodule
ifeq ($(HTS), install)
    HTS_LIB = ./htslib/libhts.a
    HTS_INCLUDE = -I./htslib
else
    # Use system-wide htslib
    HTS_LIB =
    HTS_INCLUDE ?=
    LIBS += -lhts
endif

# Default to build and link the libminimap2 submodule
ifeq ($(MINIMAP2), install)
    MINIMAP2_LIB = ./minimap2/libminimap2.a
    MINIMAP2_INCLUDE = -I./minimap2
else
    # Use system-wide htslib
    MINIMAP2_LIB =
    MINIMAP2_INCLUDE =
    LIBS += -lminimap2
endif

# Default to build and link slow5 submodule
ifeq ($(SLOW5LIB), install)
    SLOW5LIB_LIB = ./slow5lib/lib/libslow5.a
    SLOW5LIB_INCLUDE = -I./slow5lib/include/
endif

ifeq ($(ARM), 1)
    MINIMAP2_OPT += arm_neon=1
endif

# Include the src subdirectories
NP_INCLUDE = $(addprefix -I./, $(SUBDIRS))

# Add include flags
CPPFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(MINIMAP2_INCLUDE) $(NP_INCLUDE) $(EIGEN_INCLUDE) $(SLOW5LIB_INCLUDE)

# Main programs to build
PROGRAM = nanopolish
TEST_PROGRAM = nanopolish_test

.PHONY: all
all: depend $(PROGRAM)

#
# Build libhts
#
htslib/libhts.a:
	cp etc/htslib_config.h htslib/config.h
	$(MAKE) -C htslib CFLAGS="$(MAC_FLAGS)" LDFLAGS="$(MAC_FLAGS)" htslib_default_libs="-lz -lm -lbz2" NONCONFIGURE_OBJS="" || exit 255

minimap2/libminimap2.a:
	$(MAKE) -C minimap2 $(MINIMAP2_OPT) CFLAGS="$(MAC_FLAGS)" libminimap2.a || exit 255

slow5lib/lib/libslow5.a:
	$(MAKE) -C slow5lib CFLAGS="$(MAC_FLAGS)" || exit 255
#
# If this library is a dependency the user wants HDF5 to be downloaded and built.
#
lib/libhdf5.a:
	if [ ! -e hdf5-$(HDF5_VERSION).tar.gz ]; then \
		version_major_minor=`echo "$(HDF5_VERSION)" | sed -E 's/\.[0-9]+$$//'`; \
		wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$${version_major_minor}/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.gz; \
	fi

	tar -xzf hdf5-$(HDF5_VERSION).tar.gz || exit 255
	cd hdf5-$(HDF5_VERSION) && \
		./configure CFLAGS="$(MAC_FLAGS)" CPPFLAGS="$(MAC_FLAGS)" --enable-threadsafe --disable-hl --libdir=`pwd`/../lib --includedir=`pwd`/../include --prefix=`pwd`/.. || exit 255
	$(MAKE) -C hdf5-$(HDF5_VERSION) && $(MAKE) -C hdf5-$(HDF5_VERSION) install

# Download and install eigen if not already downloaded
eigen/INSTALL:
	if [ ! -e $(EIGEN_VERSION).tar.bz2 ]; then \
		wget https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.bz2; \
	fi
	tar -xjf eigen-$(EIGEN_VERSION).tar.bz2 || exit 255
	mv eigen-$(EIGEN_VERSION) eigen || exit 255

#
# Source files
#

# Find the source files by searching subdirectories
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))
EXE_SRC = src/main/nanopolish.cpp src/test/nanopolish_test.cpp

# Automatically generated object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

# Generate dependencies
.PHONY: depend
depend: .depend

.depend: $(CPP_SRC) $(C_SRC) $(EXE_SRC) $(H5_LIB) $(EIGEN_CHECK)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

# Compile objects
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS) $(H5_INCLUDE) -fPIC $<

# Link main executable
$(PROGRAM): src/main/nanopolish.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(MINIMAP2_LIB) $(H5_LIB) $(EIGEN_CHECK) $(SLOW5LIB_LIB)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(MINIMAP2_LIB) $(SLOW5LIB_LIB) $(H5_LIB) $(LIBS) $(LDFLAGS)

# Link test executable
$(TEST_PROGRAM): src/test/nanopolish_test.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(MINIMAP2_LIB) $(H5_LIB) $(SLOW5LIB_LIB)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(MINIMAP2_LIB) $(SLOW5LIB_LIB) $(H5_LIB) $(LIBS) $(LDFLAGS)

.PHONY: test
test: $(TEST_PROGRAM)
	./$(TEST_PROGRAM)

.PHONY: clean
clean:
	rm -f $(PROGRAM) $(TEST_PROGRAM) $(CPP_OBJ) $(C_OBJ) \
		src/main/nanopolish.o src/test/nanopolish_test.o
