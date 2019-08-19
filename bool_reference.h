//
// Created by sjef on 19-8-19.
//

#ifndef VPGSOLVER_BOOL_REFERENCE_H
#define VPGSOLVER_BOOL_REFERENCE_H

#include <unordered_set>
using namespace std;

struct bool_reference {
    unordered_set<int> * ref;
    int element;
    explicit bool_reference(unordered_set<int> * __ref, int element){
        ref = __ref;
        this->element = element;
    }

    explicit operator bool() const
    {
        return ref->find(element) != ref->end();
    }

    bool_reference& operator=(bool b){
        if(b)
            ref->insert(element);
        else
            ref->erase(element);
        return *this;
    }
};


#endif //VPGSOLVER_BOOL_REFERENCE_H
