CXX= mpicxx
INC=-I/usr/local/include/blaze -I/usr/local/include/eigen3 -I/software/boost
SRC=main_mpi.cpp GLP/GLP-Data/GLP_Data.cpp GLP/GLP-DisUtensils/GLP_DisUtensils.cpp GLP/GLP-Functions/GLP_functions.cpp GLP/GLP-Libraries/GLP_libraries.cpp GLP/GLP-Libraries/GLP_splines.cpp GLP/GLP-Priors/GLP_Priors.cpp GLP/GLP-Recipes/GLP_Recipes.cpp
OBJ=*.o
CXXFLAGS=-std=c++14 -O3 -DNDEBUG -mavx -pthread -fopenmp
MKLPATH = /opt/intel/oneapi/mkl/latest/lib/intel64
MKLINK= -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_gnu_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
LDFLAGS=-L/software/boost/stage/lib
LIBS= -L/opt/intel/oneapi/mkl/2022.0.2/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lgomp -lm -ldl -lmpi -lstdc++ -lboost_serialization -lboost_mpi -llapack
EXE=myfood
all: inlaplus.o inlaplus.exe

inlaplus.exe:
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS) $(LIBS) $(MKLINK)
inlaplus.o: 
	$(CXX) $(CXXFLAGS) -c $(SRC) $(INC)
clean:
	rm $(OBJ) $(EXE)

