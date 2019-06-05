cmake_minimum_required(VERSION 3.14)
project(VPGSolver)

set(CMAKE_CXX_STANDARD 14)

add_executable(VPGSolver main.cpp VPG/Vertex.cpp VPG/Vertex.h VPG/Edge.cpp VPG/Edge.h Game.cpp Game.h BDD/bddObj.h BDD/bddObj.cpp BDD/bbdToDot.cpp zlnk.cpp zlnk.h)