CXX = g++
CXXFLAGS = -g -Isrc -Itests  

SOURCES = src/*.cpp main.cpp
EXECUTABLE = ./build/main

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CXX) $(CXXFLAGS) $^ -o $@
