include ../local_config.mk
include config.mk

SHELL := /bin/bash
PATH := $(PATH):$(BASILISK)
BASILISK := $(BASILISK)
LD_LIBRARY_PATH := $(shell pwd)/$(I_Path)
export LD_LIBRARY_PATH=$(shell pwd)/$(I_Path)

RUN = mpirun -np $(PROCCESORS) ./a.out

run:
	@echo Copy paste the following in case of error:
	@echo export LD_LIBRARY_PATH=$(shell pwd)/$(I_Path)
	$(RUN) $(tmax) $(L) $(Z0) $(a) $(Re) $(ns) $(initial_level) $(max_level) $(min_level) $(save_dt) $(path)

source:
	qcc  -source -O2 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c -L$(shell pwd)/$(I_Path) -lconnector -lstdc++ -Wl,-rpath,$(shell pwd)/$(I_Path)
	mpicc -Wall -O2 -std=c99 _Circle_Ring.c -lm -lmpi -L$(shell pwd)/$(I_Path) -lconnector -lstdc++ -Wl,-rpath,$(shell pwd)/$(I_Path)

direct:
	qcc -O2 -std=c99 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c -lm -lmpi -L$(shell pwd)/$(I_Path) -lconnector -lstdc++ -Wl,-rpath,$(shell pwd)/$(I_Path)