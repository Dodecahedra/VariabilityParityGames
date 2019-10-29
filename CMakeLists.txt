cmake_minimum_required(VERSION 3.14)
project(VPGSolver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ../../bin/)
add_executable(VPGSolver_bdd main.cpp Game.cpp Game.h Algorithms/zlnkVPG.cpp Algorithms/zlnkVPG.h conf.h Algorithms/Datastructures/ConfSetExplicit.cpp Algorithms/Datastructures/ConfSetExplicit.h Algorithms/FPIte.cpp Algorithms/FPIte.h Algorithms/MBR.cpp Algorithms/MBR.h Algorithms/Datastructures/bintree.cpp Algorithms/Datastructures/bintree.h Algorithms/Datastructures/UnorderedVertexSet.cpp Algorithms/Datastructures/UnorderedVertexSet.h Algorithms/Datastructures/VectorBoolOptimized.cpp Algorithms/Datastructures/VectorBoolOptimized.h)
add_executable(VPGSolver_bdd_single main.cpp Game.cpp Game.h conf.h Algorithms/Datastructures/ConfSetExplicit.cpp Algorithms/Datastructures/ConfSetExplicit.h Algorithms/FPIte.cpp Algorithms/FPIte.h Algorithms/MBR.h Algorithms/Datastructures/UnorderedVertexSet.cpp Algorithms/Datastructures/UnorderedVertexSet.h Algorithms/Datastructures/VectorBoolOptimized.cpp Algorithms/Datastructures/VectorBoolOptimized.h Algorithms/zlnkPG.cpp Algorithms/zlnkPG.h)

add_executable(VPGSolver_explicit main.cpp Game.cpp Game.h Algorithms/zlnkVPG.cpp Algorithms/zlnkVPG.h conf.h Algorithms/Datastructures/ConfSetExplicit.cpp Algorithms/Datastructures/ConfSetExplicit.h Algorithms/FPIte.cpp Algorithms/FPIte.h Algorithms/MBR.cpp Algorithms/MBR.h Algorithms/Datastructures/UnorderedVertexSet.cpp Algorithms/Datastructures/UnorderedVertexSet.h Algorithms/Datastructures/VectorBoolOptimized.cpp Algorithms/Datastructures/VectorBoolOptimized.h)
add_executable(VPGlambdameasure main.cpp Game.cpp Game.h Algorithms/zlnkVPG.cpp Algorithms/zlnkVPG.h conf.h Algorithms/Datastructures/ConfSetExplicit.cpp Algorithms/Datastructures/ConfSetExplicit.h Algorithms/FPIte.cpp Algorithms/FPIte.h Algorithms/MBR.cpp Algorithms/MBR.h Algorithms/Datastructures/UnorderedVertexSet.cpp Algorithms/Datastructures/UnorderedVertexSet.h Algorithms/Datastructures/bool_reference.cpp Algorithms/Datastructures/bool_reference.h Algorithms/Datastructures/VectorBoolOptimized.cpp Algorithms/Datastructures/VectorBoolOptimized.h)
find_library(BDD_LIB bdd lib)

target_compile_definitions(VPGSolver_bdd PUBLIC subsetbdd)
target_link_libraries(VPGSolver_bdd LINK_PUBLIC ${BDD_LIB})

target_compile_definitions(VPGlambdameasure PUBLIC subsetbdd lambdameasure)
target_link_libraries(VPGlambdameasure LINK_PUBLIC ${BDD_LIB})

target_compile_definitions(VPGSolver_bdd_single  PUBLIC subsetbdd SINGLEMODE)
target_link_libraries(VPGSolver_bdd_single  LINK_PUBLIC ${BDD_LIB})

target_compile_definitions(VPGSolver_explicit PUBLIC subsetexplicit)
