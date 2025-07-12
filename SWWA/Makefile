CC = gcc
CFLAGS = -Wall -O3
LIBS = -lm -lgsl -lgslcblas
all: swwa_ising

swwa_ising: swwa_ising.c rng.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f swwa_ising
