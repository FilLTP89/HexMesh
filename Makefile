## -*- Makefile -*-
##
###############################################################################
#
# source files
srcfiles 	:= $(wildcard src/*.cpp) HexMesh.cpp
#
# object files
objects		:= $(patsubst %.cpp, %.o, $(srcfiles))
#################################################################################

GTS_DIR:=/usr/local
SC_DIR:=/path/to/libsc/local
GLIB_INCLUDE=-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include

# Mesquite (download from git@github.com:FilLTP89/mesquite.git)
MSQ_LIB=/path/to/mesquite_github/src
MSQ_INCLUDE:=$(addprefix -I,$(shell find /path/to/mesquite_github/src -type d -print))
H5_LIBDIR=/path/to/hdf5-seq/lib
H5_INCLUDEDIR=/path/to/hdf5-seq/include
#H5_FLAG = -L${H5_LIBDIR} ${H5_LIBDIR}/libhdf5_hl.a ${H5_LIBDIR}/libhdf5.a ${H5_LIBDIR}/libhdf5_cpp.a  -lz -lsz -ldl -lm -Wl,-rpath -Wl,${H5_LIBDIR}
H5_FLAG = -L${H5_LIBDIR} ${H5_LIBDIR}/libhdf5_hl_cpp.a ${H5_LIBDIR}/libhdf5_cpp.a ${H5_LIBDIR}/libhdf5_hl.a ${H5_LIBDIR}/libhdf5.a -lz -lsz -ldl -lm -Wl,-rpath -Wl,${H5_LIBDIR}

H5_INCLUDE=-I${H5_INCLUDEDIR}

MPI_LP :=/usr

CXX       = mpicxx.openmpi -O1 -g -std=c++11
LDFLAGS   = -L$(GTS_DIR)/lib -lgts -L$(SC_DIR)/lib -lsc -lm -lglib-2.0 $(H5_FLAG) -L$(MSQ_LIB) $(MSQ_LIB)/libmesquite.so 
CXX_FLAGS = -I$(GTS_DIR)/include -I$(SC_DIR)/include $(GLIB_INCLUDE) $(MSQ_INCLUDE) -Iinclude $(H5_INCLUDE) 

# Target: all

all: hexmesh.exe


hexmesh.exe: $(objects)
	$(CXX) $(objects) -o hexmesh.exe $(LDFLAGS)

%.o : %.cpp
	@echo "Compiling C++ "$<"..."
	$(CXX) $(CXX_FLAGS) -c $< -o $@


#### Clean target deletes all generated files ####
clean: 
	rm $(objects) hexmesh.exe


# Enable dependency checking
.KEEP_STATE:
.KEEP_STATE_FILE:.make.state.GNU-amd64-Linux

