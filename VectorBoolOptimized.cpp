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

void VectorBoolOptimized::copy_n(VectorBoolOptimized * other, int start, int n) {
    int i = start;
    int fn = (STEPSIZE - (i % STEPSIZE)) % STEPSIZE;
    if(fn > n) fn = n;
    copy_in_int((*other)[i]._M_p,(*this)[i]._M_p,i%STEPSIZE,fn);
//    for (int j = 0; j < fn ; j++)
//    {
//        if((*this)[i + j] != (*other)[i+j]){
//            copy_in_int((*other)[i]._M_p,(*this)[i]._M_p,i%STEPSIZE,fn);
//        }
//    }

    i += fn;
    if(i != start + n) {
        for (; i + STEPSIZE < start + n; i += STEPSIZE)
            *((*this)[i]._M_p) = *((*other)[i]._M_p);
        copy_in_int((*other)[i]._M_p, (*this)[i]._M_p, 0, start + n - i);
//        for (int j = 0; j < start + n - i ; j++)
//        {
//            if((*this)[i + j] != (*other)[i+j]){
//                copy_in_int((*other)[i]._M_p, (*this)[i]._M_p, 0, start + n - i);
//            }
//        }
    }
}

void VectorBoolOptimized::copy_in_int(_Bit_type *from, _Bit_type *to, int i, int n) {
    if(n == 0)
        return;
    _Bit_type mask = 1;
    if(n == STEPSIZE){
        mask = 0;
        mask -= 1;
    } else {
        mask <<= n;
        mask -= 1;
        mask <<= i;
    }
    _Bit_type imask = ~mask;
    _Bit_type v = *from;
    v |= imask;
    *to |= mask;
    *to &= v;
}

bool VectorBoolOptimized::compare_n(VectorBoolOptimized *other, int start, int n) {
    int i = start;
    int fn = (STEPSIZE - (i % STEPSIZE)) % STEPSIZE;
    if(fn > n) fn = n;
    if(!compare_in_int((*other)[i]._M_p,(*this)[i]._M_p,i%STEPSIZE,fn))
        return false;
    i += fn;
    if(i != start + n) {
        for (; i + STEPSIZE < start + n; i += STEPSIZE)
            if(*((*this)[i]._M_p) != *((*other)[i]._M_p))
                return false;
        return compare_in_int((*other)[i]._M_p, (*this)[i]._M_p, 0, start + n - i);
    }
    return true;
}

bool VectorBoolOptimized::compare_in_int(_Bit_type *from, _Bit_type *to, int i, int n) {
    if(n == 0)
        return true;
    _Bit_type mask = 1;
    mask <<= n;
    mask -= 1;
    mask <<= i;
    _Bit_type  from_masked = *from;
    _Bit_type to_masked = *to;
    from_masked &= mask;
    to_masked &= mask;
    return from_masked == to_masked;
}
