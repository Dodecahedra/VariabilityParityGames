//
// Created by sjef on 12-10-19.
//

#ifndef VPGSOLVER_VECTORBOOLOPTIMIZED_H
#define VPGSOLVER_VECTORBOOLOPTIMIZED_H

#include <vector>
#include <climits>   // CHAR_BIT
using namespace std;

#define STEPSIZE sizeof(_Bit_type) * CHAR_BIT
class VectorBoolOptimized : public vector<bool>{
public:
    using vector<bool>::vector;
    void operator&=(VectorBoolOptimized &other);
    void operator|=(VectorBoolOptimized &other);
    void operator-=(VectorBoolOptimized &other);
};


#endif //VPGSOLVER_VECTORBOOLOPTIMIZED_H
