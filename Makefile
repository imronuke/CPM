CXX=g++
CXXFLAGS=-std=c++17 -O3 -Iinclude #-D_CHECK=1 #-Wall -D_GG=1
# CXXFLAGS=-std=c++17 -Wall -O0 -Iinclude -rdynamic -fno-omit-frame-pointer -g -fsanitize=address -fsanitize=bounds -fsanitize=undefined -D_CHECK=1 #-Wall -D_GG=1
LDFLAGS=${CXXFLAGS}

# Executable target
cpm: src/mesh.o src/pugixml.o src/read_input.o src/cpm.o src/outer.o src/main.o
	${CXX} -o $@ ${LDFLAGS} $^

# Object file targets
src/mesh.o: src/mesh.cpp
src/pugixml.o: src/pugixml.cpp
src/read_input.o: src/read_input.cpp
src/cpm.o: src/cpm.cpp
src/outer.o: src/outer.cpp
src/main.o: src/main.cpp

# Suffixes and rule for compiling .cpp to .o
.SUFFIXES: .o .cpp
.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

.PHONY: clean
clean:
	rm -f src/*.o *.a *.core # clean object files, static libraries, and core dumps
