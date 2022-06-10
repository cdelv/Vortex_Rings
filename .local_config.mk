#For source code instalations do
#MY_MFEM_INSTALL_DIR = /home/wind/Documents/MFEM/mfem/build
#HYPRE_INC = -I$(MY_MFEM_INSTALL_DIR)/../../hypre/src/hypre/include
#SUNDIALS_INC = -I$(MFEM_INSTALL_DIR)/../../sundials/install/include

#For Spack instalations do
MY_MFEM_INSTALL_DIR = /opt/spack/opt/spack/linux-pop22-icelake/gcc-11.2.0/mfem-4.4.0-jttro3op37l25j7op5hhg4taglsotydj

CONFIG_MK=$(MY_MFEM_INSTALL_DIR)/share/mfem/config.mk 
GENERAL = $(MY_MFEM_INSTALL_DIR)/include/mfem/general

SHARE_DIR = NULL
PROCCESORS = 1

include $(CONFIG_MK)

