module load \
    slurm/alpine curc-quota/latest StdEnv \
    intel/2024.2.1 impi/2021.13 cmake/3.31.0 git/2.31.0 netcdf/4.9.2 openblas/0.3.28

export ESMF_ABI="64"
export ESMF_BOPT="O"
export ESMF_COMM="intelmpi"
export ESMF_COMPILER="intel"
export ESMF_DIR="/projects/yuli3660/library/esmf"
export ESMF_INSTALL_PREFIX="/projects/yuli3660/library/esmf"

export ESMF_LAPACK="openblas"
export ESMF_LAPACK_LIBPATH="/curc/sw/install/openblas/0.3.28/intel/2024.2.1/lib"
export ESMF_NETCDF="nc-config"
export ESMF_PIO="internal"
export ESMF_PNETCDF="pnetcdf-config"

export ESMFMKFILE="/projects/yuli3660/library/esmf/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk"
export ESMF_LIB="/projects/yuli3660/library/esmf/lib/libO/Linux.intel.64.intelmpi.default"
export ESMF_MOD="/projects/yuli3660/library/esmf/mod/modO/Linux.intel.64.intelmpi.default"
export ESMF_OS="Linux"

export JASPERINC="/projects/yuli3660/library/jasper/include"
export JASPERLIB="/projects/yuli3660/library/jasper/lib64"

export MITGCM_DIR="/home/yuli3660/projects/model/MITgcm"
export SKRIPS_DIR="/home/yuli3660/projects/model/scripps_kaust_model"
export WRF_DIR="/home/yuli3660/projects/model/WRF"

export SKRIPS_MPI_INC="/curc/sw/install/intel/2024.2.1/mpi/2021.13/include"
export SKRIPS_MPI_LIB="/curc/sw/install/intel/2024.2.1/mpi/2021.13/lib"
export SKRIPS_NETCDF_INCLUDE="/curc/sw/install/netcdf/4.9.2/impi/2021.13/intel/2024.2.1/include"
export SKRIPS_NETCDF_LIB="/curc/sw/install/netcdf/4.9.2/impi/2021.13/intel/2024.2.1/lib"

export LD_LIBRARY_PATH="$ESMF_LIB:$LD_LIBRARY_PATH"
