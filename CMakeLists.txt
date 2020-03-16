#-----------------------------------------------------------------------------#
#/*!
# * \file   CMakeLists.txt
# * \author Stefano Tognini
# * \brief  CMakeLists for the H5ToRoot package within the Celeritas Project.
# * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
# */
#-----------------------------------------------------------------------------#


#------------------------------------------------------------------------------
# Setting up the project
#------------------------------------------------------------------------------
project(H5ToRoot LANGUAGES CXX C)

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)


#------------------------------------------------------------------------------
# Setting environmental variables
#------------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


#------------------------------------------------------------------------------
# Finding packages
#------------------------------------------------------------------------------
find_package(ROOT REQUIRED)
find_package(HDF5 REQUIRED)
find_package(SCALE COMPONENTS Nemesis)
find_package(MPI REQUIRED COMPONENTS C)


#------------------------------------------------------------------------------
# Adding the executable, and linking it to the needed libraries
#------------------------------------------------------------------------------
add_executable(h5ToRoot main.cpp src/H5ToRoot.cpp include/H5ToRoot.hh)

target_include_directories(h5ToRoot PRIVATE
  ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(h5ToRoot PRIVATE
  ROOT::Core ROOT::Tree
  ${HDF5_LIBRARIES}
  Nemesis
  MPI::MPI_C
)


#------------------------------------------------------------------------------
# Installing the executable to the current directory
#------------------------------------------------------------------------------
install(TARGETS h5ToRoot DESTINATION .)
