//
// Created by sjef on 12-10-19.
//

#ifndef VPGSOLVER_VECTORBOOLOPTIMIZED_H
#define VPGSOLVER_VECTORBOOLOPTIMIZED_H

#include <vector>
#include <climits>   // CHAR_BIT
using namespace std;

// How many bits are stored in the underlying structure used by vector<bool>
#define STEPSIZE (sizeof(_Bit_type) * CHAR_BIT)
/**
 * vector<bool> uses integers to efficiently represent many boolean variables.
 * boolean vectors can be used to represent sets, when performing set operation (e.g. union) then we do not want to
 * iterate all bools, we simply want to the perform bitwise and operator on the underlying integers.
 *
 * Also when we want to copy or compare parts of boolean vectors we want to do so by modifying/comparing the underlying
 * integers as opposed to the bools individually
 *
 * This class performs the operations above in an efficient manner.
 *
 * The C++ boost library should have been used for this.
 */
class VectorBoolOptimized : public vector<bool>{
public:
    using vector<bool>::vector;

    void operator&=(VectorBoolOptimized &other);
    void operator|=(VectorBoolOptimized &other);
    void operator-=(VectorBoolOptimized &other);

    /**
     * Copy n boolean values from other starting at index start
     * @param other copy from
     * @param start start at
     * @param n length
     */
    void copy_n(VectorBoolOptimized * other, int start, unsigned int n);
    /**
     * Compare n boolean values from other starting at index start
     * @param other compare to
     * @param start start at
     * @param n length
     * @return true iff all boolean values from start to start+n are equal for this and other
     */
    bool compare_n(VectorBoolOptimized *  other, int start, unsigned int n);

private:
    /**
     * Copy n bits from "from" to "to" starting at the ith least significant bit
     * @param from copy from
     * @param to copy to
     * @param i start at (least significant bit)
     * @param n length
     */
    static void copy_in_int(_Bit_type * from, _Bit_type * to, int i, unsigned int n);
    /**
     * Compare n bits from "from" with "to" starting at the ith least significant bit
     * @param from copy from
     * @param to copy to
     * @param i start at (least significant bit)
     * @param n length
     * @return true iff all bits from i to i+n are equal
     */
    static bool compare_in_int(_Bit_type * from, _Bit_type * to, int i,unsigned int n);
};


#endif //VPGSOLVER_VECTORBOOLOPTIMIZED_H
