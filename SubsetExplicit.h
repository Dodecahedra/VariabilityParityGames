//
// Created by sjef on 11-6-19.
//

#ifndef VPGSOLVER_SUBSETEXPLICIT_H
#define VPGSOLVER_SUBSETEXPLICIT_H
#include <set>
#include <unordered_set>
using namespace std;

class SubsetExplicit {
public:
    SubsetExplicit();
    SubsetExplicit(int range,int bit);
    set<int> items;
    static SubsetExplicit SetFullset;
    static SubsetExplicit SetEmptyset;
    static SubsetExplicit *getFullSet();
    static SubsetExplicit *getEmptySet();
    void operator&=(SubsetExplicit &other);
    void operator|=(SubsetExplicit &other);
    void operator-=(SubsetExplicit &other);
    int operator==(SubsetExplicit &other);
    void operator=(SubsetExplicit other);
};


#endif //VPGSOLVER_SUBSETEXPLICIT_H
