
.PHONY: all

all: parma.dylib

COMPILEFLAGS := -D__PARMAPATH__=\"$(shell pwd)\"

%.o: %.cpp Makefile
	$(CXX) $(COMPILEFLAGS) -g -c $<

parma.dylib: subroutines.o # jlparma.cpp
	$(CXX) -g  -o $@ -shared -fPIC $(CXXFLAGS) $^
