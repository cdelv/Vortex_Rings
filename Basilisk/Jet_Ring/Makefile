include ../basilisk_config.mk

SHELL := /bin/bash
PATH := $(PATH):$(BASILISK)
BASILISK := $(BASILISK)


CFLAGS += -O2
RUN = mpirun -np $(PROCCESORS) ./

all:
	qcc -source -D_MPI=1 -O2 -I$(BASILISK) -L$(BASILISK) two_rings.c 


