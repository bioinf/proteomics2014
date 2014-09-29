CXX = g++
CXXFLAGS = -Wall -std=c++11

SRCDIR = src
OUTDIR = bin

EXEC = $(OUTDIR)/prg
SOURCES = $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OUTDIR)/%.o)

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXEC)

$(OBJECTS): $(OUTDIR)/%.o : $(SRCDIR)/%.cpp
	mkdir -p $(OUTDIR)
	$(CXX) $(CXXFLAGS) -c -O3 $< -o $@

clean:
	rm -rf $(OUTDIR)

.PHONY: clean
