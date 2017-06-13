NAME := SeedGenerator
SOURCES := SeedGenerator.cpp
HEADERS :=
INCLUDE := -I/usr/include/
CXXOPTIONS ?= -O6 -Wall
CXXFLAGS := -std=c++11 $(CXXOPTIONS)
LIBS := 
.PHONY: clean

$(NAME) : $(SOURCES:.cpp=.o)
	$(CXX) $^ $(LIBS) -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) $(NAME) *~
