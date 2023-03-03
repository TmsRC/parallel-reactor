SRC = src/configuration.F90 src/simulation_support.F90 src/main.F90 
LFLAGS=-lm
FFLAGS=-O3 -J build

.PHONY: archer2 local build

archer2: CC=ftn
archer2: build

cirrus: CC=gfortran
cirrus: build

build: 
	mkdir -p build
	$(CC) -o reactor $(SRC) $(FFLAGS) $(LFLAGS)
