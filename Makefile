CXX = g++
CXXFLAGS = -g -Isrc -Itests  

SOURCES = src/*.cpp src/*.hpp main.cpp
EXECUTABLE = ./build/main

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CXX) $(CXXFLAGS) $^ -o $@
