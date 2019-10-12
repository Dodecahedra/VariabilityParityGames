//
// Created by sjef on 12-10-19.
//

#include "VectorBoolOptimized.h"

void VectorBoolOptimized::operator&=(VectorBoolOptimized &other) {
    if(this->size() != other.size())
        throw "The intersection of two boolean vectors can only be taken when they have the same size";
    for(int i = 0;i<this->size();i+=STEPSIZE)
        *((*this)[i]._M_p) &= *(other[i]._M_p);
}

void VectorBoolOptimized::operator|=(VectorBoolOptimized &other) {
    if(this->size() != other.size())
        throw "The union of two boolean vectors can only be taken when they have the same size";
    for(int i = 0;i<this->size();i+=STEPSIZE)
        *((*this)[i]._M_p) |= *(other[i]._M_p);
}

void VectorBoolOptimized::operator-=(VectorBoolOptimized &other) {
    if(this->size() != other.size())
        throw "The subtraction of two boolean vectors can only be taken when they have the same size";
    for(int i = 0;i<this->size();i+=STEPSIZE)
        *((*this)[i]._M_p) &= *(other[i]._M_p) ^ *((*this)[i]._M_p);
}
