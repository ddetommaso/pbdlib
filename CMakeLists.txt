cmake_minimum_required(VERSION 2.8)
PROJECT(pbdlib)

find_package(Armadillo 3.9 REQUIRED)
set(CMAKE_MAJOR_VERSION 1)
set(CMAKE_MINOR_VERSION 0)
set(CMAKE_PATCH_VERSION 0)

include_directories(include)
include_directories(${ARMADILLO_INCLUDE_DIRS})

ADD_LIBRARY(pbd
		src/datapoints.cpp
		src/demonstration.cpp
		src/parameters.cpp
		src/mvn.cpp
		src/gmm.cpp
		src/gmr.cpp
		src/pgmm.cpp
)

################################################################################
### test
################################################################################
OPTION(PBDLIB_BUILD_TEST "Build PbDLib test" 1)
IF (PBDLIB_BUILD_TEST)
  ADD_SUBDIRECTORY(test)
ENDIF ()

################################################################################
### samples
################################################################################
OPTION(PBDLIB_BUILD_TEST "Build PbDLib samples" 1)
IF (PBDLIB_BUILD_TEST)
  ADD_SUBDIRECTORY(samples)
ENDIF ()




target_link_libraries(pbd ${ARMADILLO_LIBRARIES})
