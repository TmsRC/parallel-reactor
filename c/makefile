SRC = src/simulation_configuration.c src/main.c src/simulation_support.c src/parallel_support.c
LFLAGS=-lm
CFLAGS=-O3

.PHONY: archer2 local build

#archer2: CC=mpicc
#archer2: build

cirrus: CC=mpicc
cirrus: build

build: 
	$(CC) -o reactor $(SRC) $(CFLAGS) $(LFLAGS)

clean:
	rm ./reactor
