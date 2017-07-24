NAME := gfrdVisualizer
SOURCES := ExternalSim.cpp \
	ScreenShot.cpp \
	gfrdRenderer.cpp \
	gfrdVisualizer.cpp
HEADERS := CameraController.hpp \
	ExternalSim.hpp \
	ScreenShot.hpp \
	gfrdRenderer.hpp
INCLUDE := -I. -I../Logger -I../../GreensFunctions -I../../eGFRD -I/usr/include/
CXXOPTIONS ?= -O6 -Wall $(DEFMATSIZE)
CXXFLAGS := -std=c++14 $(CXXOPTIONS)
LIBS := -L../Logger -L../../GreensFunctions -L../../eGFRD -lrt -lglut -lGLU -lGL -lgsl -lgslcblas -leGFRD -lGreensFunctions -lLogger
.PHONY: visualizer clean

visualizer : $(NAME)

$(NAME) : $(SOURCES:.cpp=.o)
	$(CXX) $^ $(LIBS) -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) $(NAME) *~
