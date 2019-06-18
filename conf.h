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
#define Subset bdd
#define fullset bddtrue
#define emptyset bddfalse
#include "bdd.h"
#endif
#ifdef subsetexplicit
#define Subset SubsetExplicit
#define fullset SubsetExplicit::SetFullset
#define emptyset SubsetExplicit::SetEmptyset
#include "SubsetExplicit.h"
#endif

#endif //VPGSOLVER_CONF_H
