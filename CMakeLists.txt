cmake_minimum_required(VERSION 3.14)
project(VPGSolver)

set(CMAKE_CXX_STANDARD 14)

add_executable(VPGSolver main.cpp VPG/Vertex.cpp VPG/Vertex.h VPG/Edge.cpp VPG/Edge.h Game.cpp Game.h zlnk.cpp zlnk.h conf.h SubsetExplicit.cpp SubsetExplicit.h)
find_library(BDD_LIB bdd lib)
target_link_libraries(VPGSolver LINK_PUBLIC ${BDD_LIB})