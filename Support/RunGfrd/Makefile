NAME := RunGfrd
SOURCES := RunGfrd.cpp
HEADERS := RunGfrd.hpp \
	Custom.hpp \
	Simulation.hpp \
	Simulation3P.hpp \
	Equilibrium.hpp \
	mapK.hpp \
	Resume.hpp \
	Reversible.hpp \
	PowerSpectrum.hpp \
	ParseIni/SimulatorSettings.hpp \
	ParseIni/SectionBase.hpp \
	ParseIni/WorldSection.hpp \
	ParseIni/SpeciesTypeSection.hpp \
	ParseIni/SimulatorSection.hpp \
	ParseIni/ReactionRuleSection.hpp \
	ParseIni/CopyNumbersSection.hpp \
	ParseIni/ParticlesSection.hpp \
	ParseIni/ParserExceptions.hpp \
	getoptions.hpp
INCLUDE := -I. -I../Logger -I../../GreensFunctions -I../../eGFRD -I/usr/include/
CXXOPTIONS ?= -O6 -Wall $(DEFMATSIZE)
CXXFLAGS := -std=c++14 $(CXXOPTIONS)
LIBS := -L../Logger -L../../GreensFunctions -L../../eGFRD -lrt -lgsl -lgslcblas -leGFRD -lGreensFunctions -lLogger
.PHONY: samples clean

samples : $(NAME)

$(NAME) : $(SOURCES:.cpp=.o)
	$(CXX) $^ $(LIBS) -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean :
	@$(RM) $(SOURCES:.cpp=.o) $(NAME) *~
