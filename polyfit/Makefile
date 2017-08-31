INCLUDE = -Iinclude/

# compilers
cc = g++ -Wno-write-strings -O3 $(INCLUDE)

object/containers.o: src/containers.cpp include/containers.h
	$(cc) -c -o object/containers.o src/containers.cpp

polyfit: object/containers.o src/polyfit.cpp
	$(cc) -o bin/polyfit src/polyfit.cpp object/containers.o
