# life history model
EXE_LH=stress_lh.exe
CPP_LH=stress_lh.cpp

# compiler is gcc
CXX=g++
CXXFLAGS=-Wall  -O3 

all : $(EXE_LH)

$(EXE_LH) : $(CPP_LH)
	$(CXX) $(CXXFLAGS) -o $(EXE_LH) $(CPP_LH)

clean :
	rm -rf $(EXE_LH)
