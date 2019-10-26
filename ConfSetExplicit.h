//
// Created by sjef on 11-6-19.
//

#ifndef VPGSOLVER_CONFSETEXPLICIT_H
#define VPGSOLVER_CONFSETEXPLICIT_H
#include <vector>
#include "VectorBoolOptimized.h"
using namespace std;

/**
 * Decorate the VectoBoolOptimized to represent sets explicitly
 *
 * @todo: Instead of decorating, this class should extend VectorBoolOptimized.
 */
class ConfSetExplicit {
public:
    /**
     * Maximum ize of the sets
     */
    static int size;
    /**
     * Instance of SubsetExplicit containing all elements
     */
    static ConfSetExplicit SetFullset;
    /**
     * Instance of SubsetExplicit containing no elements
     */
    static ConfSetExplicit SetEmptyset;

    /**
     * Initiate an explicit set containing no items
     */
    ConfSetExplicit();
    /**
     * Initiate an explicit set according to an integer such that for every i where the ith least significant bit of the
     * integer is b we get items[i] = b
     * @param bit integer representing the which elements are in and out
     */
    explicit ConfSetExplicit(int bit);
    void operator&=(ConfSetExplicit &other);
    void operator|=(ConfSetExplicit &other);
    void operator-=(ConfSetExplicit &other);
    int operator==(ConfSetExplicit &other);
    ConfSetExplicit& operator=(ConfSetExplicit other);

    /**
     * Return the number of elements in the set
     * @return the number of elements in the set
     */
    int count();
    VectorBoolOptimized items;
};


#endif //VPGSOLVER_CONFSETEXPLICIT_H
