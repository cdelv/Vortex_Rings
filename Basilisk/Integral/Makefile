include ../local_config.mk
include config.mk

SHELL := /bin/bash
PATH := $(PATH):$(BASILISK)
BASILISK := $(BASILISK)
LD_LIBRARY_PATH := $(shell pwd)/$(I_Path)
export LD_LIBRARY_PATH=$(shell pwd)/$(I_Path)
BYellow=\033[1;33m
BCyan=\033[1;36m
BRed=\033[1;31m
NC=\033[0m

RUN = mpirun -np $(PROCCESORS) ./a.out

.PHONY: graph clean source run

all: source run

run:
	@echo -e "${BYellow}Copy paste the following in case of error:${NC}"
	@echo "export LD_LIBRARY_PATH=$(shell pwd)/$(I_Path)"
	@$(RUN) $(tmax) $(L) $(Z0) $(a) $(Re) $(ns) $(initial_level) \
	$(max_level) $(min_level) $(save_dt) $(path)

source:
	qcc  -source -O2 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c -L$(shell pwd)/$(I_Path) -lconnector -lstdc++ -Wl,-rpath,$(shell pwd)/$(I_Path)
	mpicc -Wall -O2 -std=c99 _Circle_Ring.c -lm -lmpi -L$(shell pwd)/$(I_Path) -lconnector -lstdc++ -Wl,-rpath,$(shell pwd)/$(I_Path)

direct:
	qcc -O2 -std=c99 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c -lm -lmpi -L$(shell pwd)/$(I_Path) -lconnector -lstdc++ -Wl,-rpath,$(shell pwd)/$(I_Path)

graph:
ifeq ($(SHARE_FOLDER), NULL)
	@echo -e "${BRed}No share directory.${NC}"
else
	@echo -e "${BCyan}Copying output to $(SHARED_FOLDER)${NC}"
	@rm -rf $(SHARED_FOLDER)/htg/*
	@rm -f $(SHARED_FOLDER)/output_htg.pvd
	@cp -r htg/ output_htg.pvd -t $(SHARED_FOLDER)
endif

clean:
	@echo -e "${BYellow}Deleting files ended in .out .pvd and .htg${NC}"
	@find . -type f -name "*.out" -delete
	@find . -type f -name "*.pvd" -delete
	@find . -type f -name "*.htg" -delete
