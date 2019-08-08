cmake_minimum_required(VERSION 3.14)
project(VPGSolver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ../../bin/)
add_executable(VPGSolver_bdd main.cpp Game.cpp Game.h zlnk.cpp zlnk.h conf.h SubsetExplicit.cpp SubsetExplicit.h FPIte.cpp FPIte.h MBR.cpp MBR.h Measurements/bintree.cpp Measurements/bintree.h)
add_executable(VPGSolver_bdd_single main.cpp Game.cpp Game.h zlnk.cpp zlnk.h conf.h SubsetExplicit.cpp SubsetExplicit.h FPIte.cpp FPIte.h MBR.h)

add_executable(VPGSolver_explicit main.cpp Game.cpp Game.h zlnk.cpp zlnk.h conf.h SubsetExplicit.cpp SubsetExplicit.h FPIte.cpp FPIte.h MBR.cpp MBR.h)
add_executable(VPGlambdameasure main.cpp Game.cpp Game.h zlnk.cpp zlnk.h conf.h SubsetExplicit.cpp SubsetExplicit.h FPIte.cpp FPIte.h MBR.cpp MBR.h)
find_library(BDD_LIB bdd lib)

target_compile_definitions(VPGSolver_bdd PUBLIC subsetbdd)
target_link_libraries(VPGSolver_bdd LINK_PUBLIC ${BDD_LIB})

target_compile_definitions(VPGlambdameasure PUBLIC subsetbdd lambdameasure)
target_link_libraries(VPGlambdameasure LINK_PUBLIC ${BDD_LIB})

target_compile_definitions(VPGSolver_bdd_single  PUBLIC subsetbdd SINGLEMODE)
target_link_libraries(VPGSolver_bdd_single  LINK_PUBLIC ${BDD_LIB})

target_compile_definitions(VPGSolver_explicit PUBLIC subsetexplicit)
