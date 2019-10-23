//
// Created by sjef on 12-10-19.
//

#ifndef VPGSOLVER_VECTORBOOLOPTIMIZED_H
#define VPGSOLVER_VECTORBOOLOPTIMIZED_H

#include <vector>
#include <climits>   // CHAR_BIT
using namespace std;

#define STEPSIZE (sizeof(_Bit_type) * CHAR_BIT)
class VectorBoolOptimized : public vector<bool>{
public:
    using vector<bool>::vector;
    void operator&=(VectorBoolOptimized &other);
    void operator|=(VectorBoolOptimized &other);
    void operator-=(VectorBoolOptimized &other);
    void copy_n(VectorBoolOptimized * other, int start, int n);
    bool compare_n(VectorBoolOptimized *  other, int start, int n);

    void copy_in_int(_Bit_type * from, _Bit_type * to, int i, int n);
    bool compare_in_int(_Bit_type * from, _Bit_type * to, int i, int n);
};


#endif //VPGSOLVER_VECTORBOOLOPTIMIZED_H
