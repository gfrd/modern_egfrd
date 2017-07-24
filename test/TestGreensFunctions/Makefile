NAME := TestGreensFunctions
SOURCES := tinytest.cpp \
	SphericalBesselGenerator_test.cpp \
	CylindricalBesselGenerator_test.cpp \
	Factorial_test.cpp \
	GreensFunction3D_test.cpp \
	GreensFunction3DAbs_test.cpp \
	GreensFunction3DAbsSym_test.cpp \
	GreensFunction3DRadAbs_test.cpp \
	GreensFunction3DRadInf_test.cpp \
	GslSumLevinU_test.cpp
HEADERS := ../common/tinytest.h \
	../common/tinytest_specific.cpp \
	SphericalBesselGenerator_test.hpp \
	CylindricalBesselGenerator_test.hpp \
	Factorial_test.hpp \
	GreensFunction3D_test.hpp \
	GreensFunction3DAbs_test.hpp \
	GreensFunction3DAbsSym_test.hpp \
	GreensFunction3DRadAbs_test.hpp \
	GreensFunction3DRadInf_test.hpp \
	GslSumLevinU_test.hpp
INCLUDE := -I. -I../../Support/Logger -I../../GreensFunctions -I/usr/include/
CXXOPTIONS ?= -O6 -Wall
CXXFLAGS := -std=c++14 $(CXXOPTIONS)
LIBS := -L../../Support/Logger -L../../GreensFunctions -lrt -lgsl -lgslcblas -lGreensFunctions -lLogger
.PHONY: test clean

test : $(NAME)

$(NAME) : $(SOURCES:.cpp=.o)
	$(CXX) $^ $(LIBS) -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) $(NAME) *~
