//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_ZLNKVPG_H
#define VPGSOLVER_ZLNKVPG_H


#include "../VPGame.h"

/**
 * Implementation of the recursive algorithm for VPGs.
 *
 * The algorithm allows for local solving vertex 0. The solvelocal property indicates that we can terminate early when
 * vertex 0 is won by player solvelocal. If solvelocal is -1 we are not allowed terminate early, if it is set to 2 we
 * are allowed to terminate when vertex 0 is found for either player 0 or 1.
 *
 * Sets of vertices with configurations are always represented by two variables; a vector containing a ConfSet for every
 * vertex and a vertexset such that every a vertex is in iff that vertex has a non-empty set of configurations.
 * The second variable is redundant but it allows to more conveniently check if the set of vertices is empty.
 * @todo: eveluate if using only a vector of configuration sets is more efficient.
 *
 * When solving globally solvelocal should be set to -1 initially.
 * Wehn solving locally solvelocal should be set to 2 initially.
 */
class zlnkVPG {
public:
    /**
     * Empty vertex set to check if a set of vertices is empty
     */
    static VertexSetZlnk emptyvertexset;
    /**
     * Indicate if we solve are solving local.
     *
     * If solvelocal is -1 we are not allowed terminate early, if it is set to 2 we
     * are allowed to terminate when vertex 0 is found for either player 0 or 1
     */
    int solvelocal = -1;

    /**
     * Flag to enable metrics. If set to true the algorithm will output the number of configurations we were able to
     * attract simultaneously.
     */
    static bool conf_metricoutput;
    /**
     * Measure time spent in the attractor set calculation
     */
    long attracting = 0;

    /**
     * Initialize solved
     * @param game Variability parity game
     */
    explicit zlnkVPG(VPGame * game);
    /**
     * Initialize solver for the subgame containing vertices and configuration in (bigV,vc)
     * @param game Variability parity game
     * @param bigV vertices in the game
     * @param vc configurations per vertex that are in the game
     */
    zlnkVPG(VPGame * game, VertexSetZlnk * bigV, vector<ConfSet> * vc);
    /**
     * Solve the parity game
     * Assume the winning sets are empty
     *
     * @param W0bigV Output parameter: winning set for player 0
     * @param W0vc Output parameter: winning set for player 0
     * @param W1bigV Output parameter: winning set for player 1
     * @param W1vc Output parameter: winning set for player 1
     */
    void solve(VertexSetZlnk * W0bigV, vector<ConfSet> * W0vc,VertexSetZlnk * W1bigV, vector<ConfSet> * W1vc);

protected:
    /**
     * VPG
     */
    VPGame * game;
    /**
     * Indicate what subgame we are solving
     */
    VertexSetZlnk * bigV;
    /**
     * Indicate what subgame we are solving
     */
    vector<ConfSet> * vc;


    /**
     * Attract vertices from bigA for player "player". Vertices and configurations attracted are removed from (bigV,vc),
     * so attracting creates a subgame simultaneously.
     * @param player Player for which to attract
     * @param bigA Vertices we are attracting
     * @param ac Vertices we are attracting
     */
    void attr(int player,VertexSetZlnk * bigA,vector<ConfSet> * ac);
    /**
     * Attractor implementation using a queue to maintain which vertices were attracted and for we consider predecessors
     * Runs in O(n*c^2*e) when configuration sets are explicit and O(n*c^3*e) when configuration sets are symbolic
     */
    void attrQueue(int player,VertexSetZlnk * bigA,vector<ConfSet> * ac);
    /**
     * Find the highest and lowest priority
     * @return tuple containing the lowest priority first and the highest second
     */
    tuple<int,int> getHighLowPrio();
    /**
     * Get all the vertices with a specific priority
     * @param bigA Output argument in which we place the vertices found
     * @param ac Output argument in which we place the vertices with configurations found
     * @param prio priority we are looking for
     */
    void getVCWithPrio(VertexSetZlnk *bigA, vector<ConfSet> *ac, int prio);
    /**
     * Calculate the union of (bigA,ac) and (bigB,bc) and store the result in (bigA,ac)
     * @param bigA Left hand set
     * @param ac Left hand set
     * @param bigB Right hand set
     * @param bc Right hand set
     */
    void unify(VertexSetZlnk * bigA, vector<ConfSet> *ac, VertexSetZlnk * bigB, vector<ConfSet> * bc);
    /**
     * Remove vertices in (bigA,ac) from (bigV,bc)
     * @param bigA Vertices to remove
     * @param ac Vertices to remove
     */
    void removeFromBigV(VertexSetZlnk * bigA, vector<ConfSet> *ac);
    /**
     * Remove configuration c for vertex i
     * @param i vertex from which to remove
     * @param c configurations to remove
     */
    void removeFromBigV(int i, ConfSet c);
    /**
     * Check if we are allowed to terminate early if we find vertex 0 for player "player"
     */
    bool inSolveLocal(int player);
    /**
     * Remove all configurations in C from (bigA,ac)
     * @param bigA Remove from
     * @param ac Remove from
     * @param C Set to remove
     * @return true iff (bigA,ac) is not empty
     */
    bool removeCSet(VertexSetZlnk *bigA, vector<ConfSet> *ac, ConfSet C);
};


#endif //VPGSOLVER_ZLNKVPG_H
