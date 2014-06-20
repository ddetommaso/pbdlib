cmake_minimum_required(VERSION 2.8)
PROJECT(pbdlib)

add_executable(test_datapoints test_datapoints.cpp)
target_link_libraries(test_datapoints pbd)

add_executable(test_demonstration test_demonstration.cpp)
target_link_libraries(test_demonstration pbd)

add_executable(test_demoparams test_demoparams.cpp)
target_link_libraries(test_demoparams pbd)