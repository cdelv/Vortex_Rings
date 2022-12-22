include ../local_config.mk
include config.mk

SHELL := /bin/bash
PATH := $(PATH):$(BASILISK)
BASILISK := $(BASILISK)

RUN = mpirun -np $(PROCCESORS) ./a.out

run:
	mpirun -np 1 ./a.out $(tmax) $(L) $(Z0) $(a) $(tol) $(Re) $(initial_level) $(max_level) $(threshold) $(save_dt) $(path)

source:
	qcc  -source -O2 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c
	mpicc -Wall -O2 -std=c99 _Circle_Ring.c -lm -lmpi -L$(BASILISK)/gl -I$(BASILISK)/gl -lglutils -lfb_osmesa -lGLU -lOSMesa

direct:
	qcc -O2 -std=c99 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c -lm -lmpi -L$(BASILISK)/gl -I$(BASILISK)/gl -lglutils -lfb_osmesa -lGLU -lOSMesa