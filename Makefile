CC=g++
CXXFLAGS=-I/home/pawel/opt/include/gmtl-0.6.1

test: test.o

clean:
	rm -rf *.o

.PHONY: clean
