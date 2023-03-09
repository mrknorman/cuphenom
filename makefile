#OBJS specifies which files to compile as part of the project
OBJS = ./src/cuphenom.c ./src/cuda_functions.cu

#CC specifies which compiler we're using
CC = nvcc

#COMPILER_FLAGS specifies the additional compilation options we're using

INCLUDE        = -I./include -I./io_tools/include -I./py_tools/include -I./cuda_tools/include
LAL_INCLUDE    = -I/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py39/include
PYTHON_INCLUDE = -I/home/michael.norman/.conda/envs/dragon/include/python3.10 -I/home/michael.norman/.conda/envs/dragon/lib/python3.10/site-packages/numpy/core/include/
NVIDIA_FLAGS   = -arch=native -use-fast-math -O3 -Xcompiler -march=native -Xcompiler -Ofast -forward-unknown-to-host-compiler
NVIDIA_DEBUG   = -G -Xcompiler -rdynamic

DEBUG_FLAG   = -g
WARNING_FLAG = -Wall -Wextra -Wconversion

#LINKER_FLAGS specifies the libraries we're linking against

LAL_LINK    = -L/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py39/lib/ -llal -llalsimulation
PYTHON_LINK = -L/home/michael.norman/.conda/envs/dragon/lib/python3.10/config-3.10-x86_64-linux-gnu -L/home/michael.norman/.conda/envs/dragon/lib -lpython3.10

LINKER_FLAGS = -lcrypt -lpthread -ldl -lgsl -lutil -lrt -lm -Wl,-rpath,/home/michael.norman/.conda/envs/dragon/lib -Wl,-rpath,/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py39/lib/

#OBJ_NAME specifies the name of our exectuable
OBJ_NAME   = ./bin/main.out

#This is the target that compiles our executable
all : $(OBJS)
	$(CC) $(OBJS) -pg -lcurand -g -lcufft $(INCLUDE) $(PYTHON_INCLUDE) $(LAL_INCLUDE) $(WARNING_FLAG) $(NVIDIA_DEBUG) $(NVIDIA_FLAGS) $(LAL_LINK) $(PYTHON_LINK) $(LINKER_FLAGS) -o $(OBJ_NAME)
	