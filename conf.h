//
// Created by sjef on 11-6-19.
//

#ifndef VPGSOLVER_CONF_H
#define VPGSOLVER_CONF_H

/**
 * Define subsetbdd to indicate the sets of configurations ought to be expressed as bdd
 * Define subsetexplicit to indicate the sets of configurations ought to be expressed as bdd
 * Define VertexSetZlnkIsBitVector to represent sets of vertices used by the recursive algorithm as bit vectors
 * Define VertexSetZlnkIsHashSet to represent sets of vertices used by the recursive algorithm as unordered sets
 *
 * Define VertexSetFPIte to be the type used to represent sets of vertices in the fixed-point iteration algorithm
 */
#ifdef subsetbdd
#include "bdd.h"
#define ConfSet bdd
#define fullset bddtrue
#define emptyset bddfalse
#endif
#ifdef subsetexplicit
#define ConfSet ConfSetExplicit
#define fullset ConfSetExplicit::SetFullset
#define emptyset ConfSetExplicit::SetEmptyset
#include "Algorithms/Datastructures/ConfSetExplicit.h"
#endif
#define VertexSetZlnkIsBitVector
//#define VertexSetZlnkIsHashSet
#ifdef VertexSetZlnkIsBitVector
#include "Algorithms/Datastructures/VectorBoolOptimized.h"
#define VertexSetZlnk  VectorBoolOptimized
#endif
#ifdef VertexSetZlnkIsHashSet
#include "Algorithms/Datastructures/UnorderedVertexSet.h"
#define VertexSetZlnk  UnorderedVertexSet
#endif

#define VertexSetFPIte VectorBoolOptimized
#endif //VPGSOLVER_CONF_H
