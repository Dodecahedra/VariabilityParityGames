//
// Created by sjef on 17-7-19.
//

#ifndef VPGSOLVER_BINTREE_H
#define VPGSOLVER_BINTREE_H

using namespace std;

template <class T>
class bintree {
public:
    T *value = nullptr;
    bintree<T> *left = nullptr;
    bintree<T> *right = nullptr;


    bool isLeaf();
};


#endif //VPGSOLVER_BINTREE_H
