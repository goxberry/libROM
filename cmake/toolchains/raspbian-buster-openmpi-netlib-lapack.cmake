# Toolchain assuming a Raspbian Buster (Debian 10 for Raspberry Pi)
# with all prerequisites installed via apt-get. Note that CMake
# auto-detection goes awry here -- which is a bug -- because of the
# way Debian/Raspbian installs ScaLAPACK.

# apt-get install gcc g++ gfortran cmake libblas3 liblapack3
# apt-get install libhdf5-dev libscalapack-openmpi-dev doxygen
# apt-get install libgtest-dev

# Raspbian uses a bizarre architecture - their own flavor of armhf:
# - it's a special flavor of armhf
# - normal armhf is armv7 + "hard float"
# - Raspbian's "armhf" is armv6 + "hard float"
# - "hard float" = "hardware floating point support" (VFP2)
# - contrast with "armel" = armv4 + "soft float"
# - "soft float" = "software floating point" via integer registers
# - of course, "soft float" is slower than "hard float"
# - also contrast with "aarch64" = arm64 = armv8 + "hard float"
# - aarch64 = 64-bit kernel + 64-bit userland
# - Raspbian buster uses 64-bit kernel, but 32-bit userland
set(CMAKE_LIBRARY_ARCHITECTURE arm-linux-gnueabihf)

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/bin/mpif90)

set(BLA_VENDOR Generic)
set(BLAS_ROOT /usr)
set(LAPACK_ROOT /usr)

# NOTE(oxberry1@llnl.gov, goxberry@gmail.com): A bug in
# Debian/Raspbian 10 (buster) writes the wrong installation path to
# the CMake config files installed by Netlib's ScaLAPACK. The
# (admittedly very hacky) fix is to symlink the ScaLAPACK library from
# the location where it is (correctly) installed to the path in the
# ScaLAPACK CMake config file, e.g.:
#
# sudo ln -s /usr/lib/arm-linux-gnueabihf/libscalapack-openmpi.so.2.0.2 \
#   /usr/lib/libscalapack-openmpi.so.2.0.2
set(ScaLAPACK_ROOT /usr)

set(HDF5_ROOT /usr)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Output compilation database" FORCE)
