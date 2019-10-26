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
    // for any bool a and b we should set a to 0 if b=1, we can express this by a = a & (b ^a)
    for(int i = 0;i<this->size();i+=STEPSIZE)
        *((*this)[i]._M_p) &= *(other[i]._M_p) ^ *((*this)[i]._M_p);
}

void VectorBoolOptimized::copy_n(VectorBoolOptimized * other, int start, unsigned int n) {
    int i = start;
    //start probably lays somewhere within an integer. First we copy the remaining part of that integer so we end up exactly at some index divisably by 64
    int fn = (STEPSIZE - (i % STEPSIZE)) % STEPSIZE;
    if(fn > n) fn = n;
    copy_in_int((*other)[i]._M_p,(*this)[i]._M_p,i%STEPSIZE,fn);

//debug code
//    for (int j = 0; j < fn ; j++)
//    {
//        if((*this)[i + j] != (*other)[i+j]){
//            copy_in_int((*other)[i]._M_p,(*this)[i]._M_p,i%STEPSIZE,fn);
//        }
//    }

    i += fn; // we did the first fn steps
    if(i != start + n) {
        // still work to do, copy complete integers for integers such that all boolean values represented by them are in the range of start to start+n
        for (; i + STEPSIZE < start + n; i += STEPSIZE)
            *((*this)[i]._M_p) = *((*other)[i]._M_p);
        // copy the part of the last integer
        copy_in_int((*other)[i]._M_p, (*this)[i]._M_p, 0, start + n - i);

// debug code
//        for (int j = 0; j < start + n - i ; j++)
//        {
//            if((*this)[i + j] != (*other)[i+j]){
//                copy_in_int((*other)[i]._M_p, (*this)[i]._M_p, 0, start + n - i);
//            }
//        }
    }
}

void VectorBoolOptimized::copy_in_int(_Bit_type *from, _Bit_type *to, int i, unsigned int n) {
    if(n == 0)
        return;
    // Create mask such that from the right to the left it contains i 0s, n 1s and then only 0s
    // The bits that are 1 in the mask are the bits that need to be copied
    _Bit_type mask = 1;
    if(n == STEPSIZE){//if n is STEPSIZE than we need mask to be completely 1s
        mask = 0;
        mask -= 1;
    } else {
        mask <<= n;
        mask -= 1;
        mask <<= i;
    }
    _Bit_type imask = ~mask;
    // create v such v is equal to from for all bits that need to be copied and 1 for all other bits
    _Bit_type v = *from;
    v |= imask;

    // Set all bits that need to be copied to 1 in to
    *to |= mask;
    // Copy all bits that need to be copied from v to to
    *to &= v;
}

bool VectorBoolOptimized::compare_n(VectorBoolOptimized *other, int start,unsigned int n) {
    // similar to copy_n


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

bool VectorBoolOptimized::compare_in_int(_Bit_type *from, _Bit_type *to, int i, unsigned int n) {
    // similar to copy_in_int


    if(n == 0)
        return true;
    _Bit_type mask = 1;
    if(n == STEPSIZE){//if n is STEPSIZE than we need mask to be completely 1s
        mask = 0;
        mask -= 1;
    } else {
        mask <<= n;
        mask -= 1;
        mask <<= i;
    }
    _Bit_type  from_masked = *from;
    _Bit_type to_masked = *to;
    from_masked &= mask;
    to_masked &= mask;
    return from_masked == to_masked;
}
