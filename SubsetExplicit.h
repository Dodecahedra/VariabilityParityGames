//
// Created by sjef on 11-6-19.
//

#ifndef VPGSOLVER_SUBSETEXPLICIT_H
#define VPGSOLVER_SUBSETEXPLICIT_H
#include <vector>
#include "VectorBoolOptimized.h"
using namespace std;

class SubsetExplicit {
public:
    SubsetExplicit();
    SubsetExplicit(int bit);
    VectorBoolOptimized items;
    static SubsetExplicit SetFullset;
    static SubsetExplicit SetEmptyset;
    static SubsetExplicit *getFullSet();
    static SubsetExplicit *getEmptySet();
    static int size;
    void operator&=(SubsetExplicit &other);
    void operator|=(SubsetExplicit &other);
    void operator-=(SubsetExplicit &other);
    int operator==(SubsetExplicit &other);
    void operator=(SubsetExplicit other);
    int count();
};


#endif //VPGSOLVER_SUBSETEXPLICIT_H
