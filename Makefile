include local_config.mk

#Compiling parameters
#path en spack
SANITFLAGS= -fsanitize=address -fsanitize=leak -fsanitize=undefined
VALGRINDFLAGS= --tool=memcheck --track-origins=yes --leak-check=full -s
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS) -I$(GENERAL) $(BOOST) $(HYPRE_INC) $(SUNDIALS_INC) #-g #$(SANITFLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard Navier/*.cpp)
DEPENDENCIES = $(SOURCES:Navier/%.cpp=.objects/%.o)

.PHONY: all main mesh graph clean oclean
.PRECIOUS: %.x $(DEPENDENCIES)

all: ellipse

test:
	@echo $(MFEM_INSTALL_DIR)
	@echo $(GENERAL)

%: %.x
	@echo -e 'Running program ... \n'
	@$(RUN)$< #$(VALGRINDFLAGS)

# ellipse: ellipse.x
# 	@echo -e 'Running program ... \n'
# 	@$(RUN)$< #$(VALGRINDFLAGS)

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/vortex_rings
	@cp -r results $(SHARE_DIR)/vortex_rings
	@echo -e 'Done!'
endif

%.x: $(DEPENDENCIES) Code/%.cpp
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

mesh: settings/struct_cube.geo
	@gmsh $< -format msh2 -o results/mesh.msh -3 > /dev/null
	@echo 'Mesh created!'

.objects/%.o: Navier/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

clean:
	@rm -rf results/graph/ *.x

oclean:
	@rm -rf .objects/*.o
