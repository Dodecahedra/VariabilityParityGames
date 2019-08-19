//
// Created by sjef on 19-8-19.
//

#ifndef VPGSOLVER_UNORDEREDVERTEXSET_H
#define VPGSOLVER_UNORDEREDVERTEXSET_H

#include <unordered_set>
#include "bool_reference.h"
using namespace std;

class UnorderedVertexSet : public unordered_set<int> {
public:
    bool_reference
    operator[](int __n)
    {
        return bool_reference(this, __n);
    }
    UnorderedVertexSet();
    UnorderedVertexSet(int size);
};


#endif //VPGSOLVER_UNORDEREDVERTEXSET_H
