CXX=g++
CXXFLAGS=-std=c++17 -O3 -Iinclude #-D_CHECK=1 #-Wall -D_GG=1
# CXXFLAGS=-std=c++17 -Wall -O0 -Iinclude -rdynamic -fno-omit-frame-pointer -g -fsanitize=address -fsanitize=bounds -fsanitize=undefined -D_CHECK=1 #-Wall -D_GG=1
LDFLAGS=${CXXFLAGS}

# Executable target
cpm: mesh.o pugixml.o read_input.o cpm.o outer.o main.o
	${CXX} -o $@ ${LDFLAGS} $^

# Object file targets
mesh.o: mesh.cpp
pugixml.o: pugixml.cpp
read_input.o: read_input.cpp
cpm.o: cpm.cpp
outer.o: outer.cpp
main.o: main.cpp

# Suffixes and rule for compiling .cpp to .o
.SUFFIXES: .o .cpp
.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.a *.core # clean object files, static libraries, and core dumps
