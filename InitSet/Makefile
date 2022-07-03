include ../local_config.mk

#Compiling parameters
#path en spack
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS) -I$(GENERAL) $(HYPRE_INC) $(SUNDIALS_INC)
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

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/elasticidad
	@cp -r results $(SHARE_DIR)/elasticidad
	@echo -e 'Done!'
endif

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

mesh: current_mesh.geo
	gmsh $< -format msh2 -o mesh.msh -3 > /dev/null

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

clean:
	rm -r results/graph/

oclean:
	@rm -rf .objects/*.o
