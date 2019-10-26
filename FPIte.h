//
// Created by sjef on 11-7-19.
//

#ifndef VPGSOLVER_FPITE_H
#define VPGSOLVER_FPITE_H

#include "conf.h"
#include "Game.h"
#include <vector>
using namespace std;

/**
 * Implementation of the fixed-point iteration algorithm for parity games.
 *
 * The algorithm is modified to work more efficiently when for some vertices it is already known who the winner is. The
 * algorithm takes P0 and VP1 such that vertices in P0 are known to be won by player 0 and vertices not in VP1 are known
 * to be own by player 1. VP1 stands for V without P1. The algorithm calculates the winning set for player 0 (W0), we
 * know that P0 <= W0 <= VP1 (where <= is the standard subset or equal operation).
 *
 * The fixed-point variables are initialized such that least fixed-point variables start at P0 and greatest fixed-point
 * variables start at VP1. Furthermore we do not consider vertices in P0 or outside VP1 in the diamond box operation.
 *
 * When P0 = emptyset and VP1 = V then the algorithm behaves identical to the original fixed-point algorithm.
 *
 * When P0 and VP1 are incorrect then no guarantees are made by the algorithm.
 *
 * Finally, the algorithm assumes that vertices are sorted by parity first and priority second.
 */
class FPIte {
public:
    /**
     * Indicates if the algorithm can terminate early when the local vertex is found
     */
    bool solvelocal = false;
    /**
     * Measure how many vertices are considered. Only used for status reporting
     */
    long verticesconsidered = 0;
    /**
     * Measure how many diamond box operators are executed. Only used for status reporting
     */
    long dbs_executed = 0;
    /**
     * Winning set for player 0
     */
    VertexSetFPIte * W0;


    /**
     * Initialize the algorithm with no vertices that are known the be won by player 0 or 1
     * @param Parity game where vertices are sorted by parity and priority
     */
    explicit FPIte(Game * game);
    /**
     * Initialize the algoritm
     * @param game Parity game where vertices are sorted by parity and priority
     * @param P0 Set of vertices that are known to be won by player 0
     * @param VP1 Set of vertices, such that any vertex not in VP1 is known to be won by player 1
     * @param W0 Set of vertices in which the winning set will be placed
     */
    FPIte(Game * game, VertexSetFPIte * P0, VertexSetFPIte * VP1, VertexSetFPIte * W0);
    ~FPIte();

    /**
     * Calculate the winning set for player 0 and store in in W0
     */
    void solve();

    /**
     * Set P0
     * @param P0string string of integers indicating the vertex numbers, separated by ,
     */
    void setP0(char * P0string);
    /**
     * Set VP1 to P0, indicating that P0 is the exact winning set of the game.
     * This function is to test the performance increase when using P0, it is not integral to the algorithm itself.
     */
    void P0IsFull();
    /**
     * Randomly choose n distinct vertices and remove them from P0 when in P0 or add them to VP1 when not in VP1.
     * Effectively decreasing the number of vertices from which the winner is known by n.
     *
     * This function is to test the performance increase when using P0, it is not integral to the algorithm itself.
     * @param n number of vertices to remove from the already solved set
     */
    void unassist(int n);

protected:
    /**
     * Parity game where vertices are sorted by parity and priority.
     */
    Game * game;
    /**
     * Set of vertices that are known to be won by player 0
     */
    VertexSetFPIte * P0;
    /**
     * Set of vertices, such that any vertex not in VP1 is known to be won by player 1
     */
    VertexSetFPIte * VP1;
    /**
     * Number of distinct priorities
     */
    int d;
    /**
     * For every vertex we maintain the number of successors that are in ZZ for vertices with owner 0 and the number of
     * successor that are not in ZZ for vertices with owner 1.
     */
    vector<int> edgecount;


    /**
     * Iteration variables. We represent all variables using a single VertexSet. Iteration variable Z_i is only ever
     * used to check if vertices with priority i are in. So for Z_i we only have to store vertices with priority i,
     * because evey vertex has a single priority we can simply maintain a single vertex set presenting all Z_i's.
     */
    VertexSetFPIte ZZ;

    /**
     * Indicates if a vertex was in ZZ in the previous diamondbox operation
     */
    VertexSetFPIte targetWasIn;


    /**
     * Initialize all iteration variables with parity (i % 2) and i >= priority >= ie
     * @param i Initialize iteration variables from
     * @param ie Initialize iteration variables to
     */
    void init(int i, int ie);
    /**
     * DiamondBox operator that initializes the edge count variable.
     * @param Z Result
     */
    void diamondbox(VertexSetFPIte * Z);
    /**
     * DiamondBox operator that uses and maintain the edge count variable.
     * Only the predecessors of vertices with priority <= maxprio are considered.
     * @param Z Results
     * @param maxprio Maximum priority of vertices for which we consider the predecessors
     */
    void diamondbox(VertexSetFPIte * Z, int maxprio);
    /**
     * Copy all vertices with priority p from ZP to Z
     * @param Z copy to
     * @param ZP copy from
     * @param p priority of vertices to copy
     */
    void copyWithPrio(VertexSetFPIte * Z, VertexSetFPIte * ZP, int p);
    /**
     * Copy all vertices with sp <= priority <= ep and parity (sp % 2) from ZP to Z
     * @param Z copy to
     * @param ZP copy from
     * @param sp start priority of vertices to copy
     * @param ep end priority of vertices to copy
     */
    void copyWithPrio(VertexSetFPIte * Z, VertexSetFPIte * ZP, int sp, int ep);
    /**
     * Compare vertices values of Z and ZP for vertices with priority p
     * @param Z compare 1
     * @param ZP compare 2
     * @param p priority of vertices to compare
     * @return True iff any vertex v with priority p is in Z iff v is in ZP
     */
    bool compareWithPrio(VertexSetFPIte * Z, VertexSetFPIte * ZP, int p);
};


#endif //VPGSOLVER_FPITE_H
