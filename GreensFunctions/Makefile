NAME := GreensFunctions
SOURCES := CylindricalBesselGenerator.cpp \
	freeFunctions.cpp \
	GreensFunction1DAbsAbs.cpp \
	GreensFunction1DAbsSinkAbs.cpp \
	GreensFunction1DRadAbs.cpp \
	GreensFunction2DAbsSym.cpp \
	GreensFunction2DRadAbs.cpp \
	GreensFunction3D.cpp \
	GreensFunction3DAbs.cpp \
	GreensFunction3DAbsSym.cpp \
	GreensFunction3DRadAbs.cpp \
	GreensFunction3DRadInf.cpp \
	GreensFunction3DSym.cpp \
	helperFunctionsGf.cpp \
	SphericalBesselGenerator.cpp
HEADERS := DefsGf.hpp GreensFunction.hpp PairGreensFunction.hpp $(SOURCES:.cpp=.hpp)
INCLUDE := -I../Support/Logger -I/usr/include/
CXXOPTIONS ?= -O6 -Wall
CXXFLAGS := -std=c++14 -fPIC $(CXXOPTIONS)
.PHONY: lib clean

lib : lib$(NAME).so

lib$(NAME).so : $(SOURCES:.cpp=.o) 
	$(CXX) -shared -Wl,-soname,lib$(NAME).so $? -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) lib$(NAME).so *~
