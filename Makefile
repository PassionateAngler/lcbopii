CC=g++
CXXFLAGS=-I/home/pawel/opt/include/gmtl-0.6.1

test: test.o atom.o lcbopii.o
#atom.o: atom.cpp

all: test

clean:
	rm -rf *.o test atom.o lcbopii.o

.PHONY: clean
