
.PHONY: all

all: parma.dylib

%.o: %.cpp
	$(CXX) -g -c $<

parma.dylib: subroutines.o # jlparma.cpp
	$(CXX) -g  -o $@ -shared -fPIC $(CXXFLAGS) $^
