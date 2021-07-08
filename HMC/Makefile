# Basic makefile
CC = gcc
OPT = -Wall -O3 -std=c99

LIBC = -lm

PROGS = phi4

all: $(PROGS)

DEPS = phi4.h ranlxd.h

OBJ = ranlxd.o utils.o hmc.o measure.o phi4.o

.c.o: phi4.h
	$(CC) $(OPT) -c $*.c

phi4: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(LIBC)

.PHONY: clean

clean:
	rm -f $(PROGS) $(OBJ)
