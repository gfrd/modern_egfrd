NAME := genBesselTables
SOURCES := genBesselTables.cpp
HEADERS :=
INCLUDE := -I/usr/include/
CXXOPTIONS ?= -O6 -Wall
CXXFLAGS := -std=c++11 $(CXXOPTIONS)
LIBS := -lgsl -lgslcblas
.PHONY: tables clean

tables: CylindricalBesselTable.bin SphericalBesselTable.bin

$(NAME) : $(SOURCES:.cpp=.o)
	$(CXX) $^ $(LIBS) -o $@

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	@$(RM) $(SOURCES:.cpp=.o) $(NAME) *~ *BesselTable.bin

SphericalBesselTable.bin : $(NAME)
	@./$(NAME) -s $@

CylindricalBesselTable.bin : $(NAME)
	@./$(NAME) -c $@

