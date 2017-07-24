NAME := TestGFRD
SOURCES := tinytest.cpp \
	BDSimulator_test.cpp \
	Box_test.cpp \
	Cylinder_test.cpp \
	Disk_test.cpp \
	EventScheduler_test.cpp \
	GenericIterator_test.cpp \
	Identifier_test.cpp \
	MatrixSpace_test.cpp \
	Model_test.cpp \
	ParticleContainer_test.cpp \
	ParticleSimulator_test.cpp \
	Plane_test.cpp \
	Persistence_test.cpp \
	Single_test.cpp \
	Sphere_test.cpp \
	TestSimple.cpp \
	Vector2_test.cpp \
	Vector3_test.cpp \
	World_test.cpp
HEADERS := ../common/tinytest.h \
	../common/tinytest_specific.cpp \
	BDSimulator_test.hpp \
	Box_test.hpp \
	Cylinder_test.hpp \
	Disk_test.hpp \
	EventScheduler_test.hpp \
	GenericIterator_test.hpp \
	Identifier_test.hpp \
	MatrixSpace_test.hpp \
	Model_test.hpp \
	ParticleContainer_test.hpp \
	ParticleSimulator_test.hpp \
	Plane_test.hpp \
	Persistence_test.hpp \
	Single_test.hpp \
	Sphere_test.hpp \
	TestSimple.hpp \
	Vector2_test.hpp \
	Vector3_test.hpp \
	World_test.hpp
INCLUDE := -I. -I../../Support/Logger -I../../GreensFunctions -I../../eGFRD -I/usr/include/
CXXOPTIONS ?= -O6 -Wall $(DEFMATSIZE) 
CXXFLAGS := -std=c++14 $(CXXOPTIONS)
LIBS := -L../../Support/Logger -L../../GreensFunctions -L../../eGFRD -lrt -lgsl -lgslcblas -leGFRD -lGreensFunctions -lLogger
.PHONY: test clean

test : $(NAME)

$(NAME) : $(SOURCES:.cpp=.o)
	$(CXX) $^ $(LIBS) -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) $(NAME) *~
