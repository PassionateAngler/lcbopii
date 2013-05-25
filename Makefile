CC=g++
CXXFLAGS=-std=c++0x -I/home/pawel/opt/include/gmtl-0.6.1

test:  lcbopii.o atom.o test.o printers.o
#atom.o: atom.cpp
#lcbopii.o: lcbopii.cpp

all: test

clean:
	rm -rf *.o test atom.o lcbopii.o

.PHONY: clean
