#Computer parameters
MFEM_INSTALL_DIR = /opt/spack/opt/spack/linux-pop21-icelake/gcc-11.2.0/mfem-4.3.0-wor2wndreilv3jrwzqxm5a54cvqgrvaq
SHARE_DIR = NULL
PROCCESORS = 1

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
