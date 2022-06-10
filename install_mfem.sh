#Exit when any command fails
set -e
#Print commands 
set -x 

#Directories PATHS
INSTALL_DIR=$(pwd)
SUNDIALS_DIR=${INSTALL_DIR}/sundials
HYPRE_DIR=${INSTALL_DIR}/hypre
METIS_DIR=${INSTALL_DIR}/metis
MFEM_INSTALL_DIR=${INSTALL_DIR}/mfem/build

#Get Hypre source code
wget -c https://github.com/hypre-space/hypre/archive/refs/tags/v2.24.0.tar.gz -O hypre.tar.gz
tar -xf hypre.tar.gz
mv hypre-2.24.0 hypre

#Get Metis 5 source code
wget -c http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz -O metis.tar.gz
tar -xf metis.tar.gz
mv metis-5.1.0 metis

#Get SUNDIALS source code
wget -c https://github.com/LLNL/sundials/releases/download/v5.7.0/sundials-5.7.0.tar.gz -O sundials.tar.gz
tar -xf sundials.tar.gz
mv sundials-5.7.0 sundials

#Get MFEM source code
git clone https://github.com/mfem/mfem.git

######################
### INSTALL HYPRE
######################
cd ${HYPRE_DIR}/src
./configure --disable-fortran
make -j 8
make install

######################
### INSTALL METIS
######################
cd ${METIS_DIR}
make config
make -j 8
mkdir lib
ln -s ../build/Linux-x86_64/libmetis/libmetis.a lib

######################
### INSTALL SUNDAILS
######################
cd ${SUNDIALS_DIR}
mkdir build
cd ${SUNDIALS_DIR}/build
cmake -DCMAKE_INSTALL_PREFIX=${SUNDIALS_DIR}/install -DEXAMPLES_INSTALL_PATH=${SUNDIALS_DIR}/install/examples -DENABLE_MPI:BOOL=ON ${SUNDIALS_DIR}
make install -j 8

######################
### INSTALL MFEM
######################
mkdir ${MFEM_INSTALL_DIR}
cd ${MFEM_INSTALL_DIR}
cmake -DMFEM_USE_MPI:BOOL=ON -DMFEM_USE_METIS:BOOL=ON -DMFEM_ENABLE_MINIAPPS:BOOL=ON -DMFEM_USE_SUNDIALS:BOOL=ON -DMFEM_USE_ZLIB:BOOL=ON -DHYPRE_DIR=${HYPRE_DIR}/src/hypre -DMETIS_DIR=${METIS_DIR} -DSUNDIALS_DIR=${SUNDIALS_DIR}/install -DCMAKE_INSTALL_NAME_DIR=${MFEM_INSTALL_DIR}/lib -DCMAKE_INSTALL_PREFIX=${MFEM_INSTALL_DIR} -DCMAKE_INSTALL_RPATH=${MFEM_INSTALL_DIR}/lib ${INSTALL_DIR}/mfem

make -j 8
make examples -j 8
make miniapps -j 8
make tests -j 8
make install
make test





