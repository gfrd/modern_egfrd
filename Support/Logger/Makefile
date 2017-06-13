NAME := Logger
SOURCES := Logger.cpp
HEADERS := $(SOURCES:.cpp=.hpp) bitmask_operators.hpp
INCLUDE := 
CXXOPTIONS ?= -O6 -Wall
CXXFLAGS := -std=c++14 -fPIC $(CXXOPTIONS)
.PHONY: lib clean

lib: lib$(NAME).so

lib$(NAME).so : $(SOURCES:.cpp=.o) 
	$(CXX) -shared -Wl,-soname,lib$(NAME).so $^ -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) lib$(NAME).so *~
