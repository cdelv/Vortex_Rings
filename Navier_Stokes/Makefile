include ../local_config.mk

#Compiling parameters
#path en spack
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS) -I$(GENERAL)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main mesh graph clean oclean

all: main

test:
	@echo $(MFEM_INSTALL_DIR)
	@echo $(GENERAL)

main: main.x
	@echo -e 'Running program ... \n'
	@$(RUN)$<

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

mesh: settings/cube.geo
	gmsh settings/cube.geo -format msh2 -o mesh.msh -3 > /dev/null

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

oclean:
	@rm -rf .objects/*.o
