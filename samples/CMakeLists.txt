cmake_minimum_required(VERSION 2.8)
PROJECT(pbdlib)

add_executable(learn_gmm learn_gmm.cpp)
target_link_libraries(learn_gmm pbd)

add_executable(repro_gmr repro_gmr.cpp)
target_link_libraries(repro_gmr pbd)

add_executable(learn_pgmm learn_pgmm.cpp)
target_link_libraries(learn_pgmm pbd)
