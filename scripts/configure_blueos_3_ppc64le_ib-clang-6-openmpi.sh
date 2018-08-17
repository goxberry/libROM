#!/usr/bin/env bash

#NOTE: Spectrum MPI is IBM's prepackaged version of OpenMPI
#NOTE: Extra LD flags are what is used by /usr/tcetmp/bin/xlflang -v
#      when called on a Fortran 90 source file. ALCF recommends
#      including these flags when symbols like _xldipow are not found
#      See https://wiki.alcf.anl.gov/old/index.php/Compiling_and_linking

export PKG_ROOT=/usr/gapps/bdiv/blueos_3_ppc64le_ib/clang-6-openmpi
export SPECTRUM_PATH=/usr/tce/packages/spectrum-mpi
export SPECTRUM_CLANG6=${SPECTRUM_PATH}/spectrum-mpi-rolling-release-clang-6.0.0
export COMPILER_ROOT=${SPECTRUM_CLANG6}/bin
export LAPACK_LDFLAGS="-L${PKG_ROOT}/lapack/3.5.0/lib -llapack -lblas"

./configure \
    --with-CXX=${COMPILER_ROOT}/mpicxx \
    --with-FC="/usr/tcetmp/bin/xlflang" \
    --with-lapack=${PKG_ROOT}/lapack/3.5.0 \
    --with-lapack-libs="${LAPACK_LDFLAGS}" \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=${PKG_ROOT}/gtest/1.8.0 \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"

# FC="/usr/tcetmp/bin/xlflang" ./configure --with-CXX="/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-6.0.0/bin/mpicxx" --with-lapack-libs="-L/usr/gapps/bdiv/blueos_3_ppc64le_ib/clang-6-openmpi/lapack/3.5.0/lib -llapack -lblas -L/usr/tce/packages/clang/xlflang-coral-2017.06.27/ibm/xlf/lib -lxl" --with-hdf5=/usr/gapps/bdiv/blueos_3_ppc64le_ib/clang-6-openmpi/hdf5/1.8.10p1 --with-zlib=/usr/gapps/bdiv/blueos_3_ppc64le_ib/clang-6-openmpi/zlib/1.2.3/lib --with-gtest=/usr/gapps/bdiv/blueos_3_ppc64le_ib/clang-6-openmpi/gtest/1.8.0 --with-elemental=no --enable-opt=yes --enable-debug=no --with-extra-ld-flags="-L/usr/tce/packages/clang/xlflang-coral-2017.06.27/ibm/xlf/lib -lxlf90_r -lxlfmath -lm -lxlmbif -lxlopt -lxl"  "$@"

