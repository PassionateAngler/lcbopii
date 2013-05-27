CC=g++
CXXFLAGS=-std=c++0x -I/usr/include/eigen3
#CXXFLAGS=-std=c++0x 
#LDLIBS=-leigen

test:  lcbopii.o atom.o test.o printers.o
#atom.o: atom.cpp
#lcbopii.o: lcbopii.cpp

all: test

clean:
	rm -rf *.o test atom.o lcbopii.o

.PHONY: clean
