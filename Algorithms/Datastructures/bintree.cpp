//
// Created by sjef on 17-7-19.
//

#include "bintree.h"

template<class T>
bool bintree<T>::isLeaf() {
    return left == nullptr && right == nullptr ;
}
