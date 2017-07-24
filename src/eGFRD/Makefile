NAME := eGFRD
SOURCES := \
	EGFRDSimulator.cpp \
	Matrix4.cpp \
	Multi.cpp \
	Pair.cpp \
	PlanarSurface.cpp \
	Persistence.cpp \
	Single.cpp \
	Vector2.cpp \
	Vector3.cpp
HEADERS := BDPropagator.hpp \
	BDSimulator.hpp \
	Box.hpp \
	ConnectivityContainer.hpp \
	convert.hpp \
	CuboidalRegion.hpp \
	Cylinder.hpp \
	CylindricalSurface.hpp \
	DefsEgfrd.hpp \
	Disk.hpp \
	DiskSurface.hpp \
	Domain.hpp \
	DomainID.hpp \
	Event.hpp \
	EventID.hpp \
	EventScheduler.hpp \
	exceptions.hpp \
	genericIterator.hpp \
	GreenFunctionHelpers.hpp \
	helperFunctions.hpp \
	Identifier.hpp \
	InteractionRule.hpp \
	makeString.hpp \
	MatrixSpace.hpp \
	Model.hpp \
	Particle.hpp \
	ParticleContainer.hpp \
	ParticleContainerImpl.hpp \
	ParticleID.hpp \
	ParticleSimulator.hpp \
	Plane.hpp \
	randomNumberGenerator.hpp \
	Reactants.hpp \
	ReactionRecorder.hpp \
	ReactionRule.hpp \
	ReactionRuleCollection.hpp \
	ReactionRuleID.hpp \
	Region.hpp \
	SerialIDGenerator.hpp \
	Shell.hpp \
	ShellCreateUtils.hpp \
	ShellID.hpp \
	SpeciesType.hpp \
	SpeciesTypeID.hpp \
	Sphere.hpp \
	SphericalSurface.hpp \
	Structure.hpp \
	StructureContainer.hpp \
	StructureID.hpp \
	StructureType.hpp \
	StructureTypeID.hpp \
	Surface.hpp \
	Transaction.hpp \
	World.hpp \
	$(SOURCES:.cpp=.hpp)
INCLUDE := -I../GreensFunctions -I../Support/Logger -I/usr/include/
CXXOPTIONS ?= -O6 -Wall $(DEFMATSIZE)
CXXFLAGS := -std=c++14 -fPIC $(CXXOPTIONS)
LIBS := -L../../Support/Logger -L../../GreensFunctions -lgsl -lgslcblas -lGreensFunctions -lLogger
.PHONY: lib clean

lib : lib$(NAME).so

lib$(NAME).so : $(SOURCES:.cpp=.o) 
	$(CXX) -shared -Wl,-soname,lib$(NAME).so $? -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) lib$(NAME).so *~
