//
// Created by sjef on 11-6-19.
//

#ifndef VPGSOLVER_CONF_H
#define VPGSOLVER_CONF_H

//#define SINGLEMODE
//#define randombddorder

//#define subsetbdd
//#define subsetexplicit

#ifdef subsetbdd
#include "bdd.h"
#define Subset bdd
#define fullset bddtrue
#define emptyset bddfalse
#endif
#ifdef subsetexplicit
#define Subset SubsetExplicit
#define fullset SubsetExplicit::SetFullset
#define emptyset SubsetExplicit::SetEmptyset
#include "SubsetExplicit.h"
#endif
#define VertexSetZlnkIsBitVector
//#define VertexSetZlnkIsHashSet
#ifdef VertexSetZlnkIsBitVector
#include "VectorBoolOptimized.h"
#define VertexSetZlnk  VectorBoolOptimized
#endif
#ifdef VertexSetZlnkIsHashSet
#include "UnorderedVertexSet.h"
#define VertexSetZlnk  UnorderedVertexSet
#endif
#define VertexSet vector<bool>
#endif //VPGSOLVER_CONF_H
