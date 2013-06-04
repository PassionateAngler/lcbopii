CC=g++
CXXFLAGS=-g -std=c++0x -I/usr/include/eigen3 -Wall # -D_GLIBCXX_DEBUG 
#CXXFLAGS=-std=c++0x 
#LDLIBS=-leigen

#test:  lcbopii.o atom.o test.o printers.o
test_tabII: lcbopii.o atom.o test_tabII.o
atom.o: atom.cpp
lcbopii.o: lcbopii.cpp

all: test_tabII

clean:
	rm -rf *.o test atom.o lcbopii.o

.PHONY: clean
