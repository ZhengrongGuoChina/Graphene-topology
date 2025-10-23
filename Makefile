# Graphene-topology: A topology-based vector library for generating graphene-derived nanostructures
#
# Author   : Zhengrong Guo (Yulin University)
# Email    : zhengrong_guo@yulinu.edu.cn
# Copyright: 2023-2025 Zhengrong Guo

Cxx = g++ #clang++

INC = -I./include
VPATH = ./include

CFLAGS = $(INC) -std=c++20 -O3

LDLIBS += #-L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib

HeadFiles = grapology_error.hpp grapology_set.hpp grapology_vector3.hpp grapology_vertex.hpp grapology_vector.hpp grapology_pool.hpp grapology_locate.hpp grapology_topology.hpp grapology_cap.hpp grapology.h 
CompFiles = src/main.cpp

grapology: $(HeadFiles) $(CompFiles)
	$(Cxx) $(CFLAGS) $(LDLIBS) $(CompFiles) -o grapology

clean:
	rm -f grapology
