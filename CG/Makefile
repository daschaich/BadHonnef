# Basic makefile with LAPACK linking
SHELL=/bin/sh
.SUFFIXES: .o .f90 .f

CC = gcc
OPT = -Wall -O3 -std=c99 -fcx-fortran-rules

LIBC = -llapack -lblas -lm

PROGS = graphene

all: $(PROGS)

DEPS = graphene.h ranlxd.h

OBJ = ranlxd.o utils.o hmc.o graphene.o

.c.o: graphene.h
	$(CC) $(OPT) -c $*.c

graphene: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(LIBC)

.PHONY: clean

clean:
	rm -f $(PROGS) $(OBJ)
