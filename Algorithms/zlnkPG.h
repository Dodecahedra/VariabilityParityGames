//
// Created by sjef on 26-10-19.
//

#ifndef VPGSOLVER_ZLNKPG_H
#define VPGSOLVER_ZLNKPG_H


#include "Game.h"

/**
 * Implementation of Zielonka's recursive algorithm.
 *
 * The algorithm allows for local solving vertex 0. The solvelocal property indicates that we can terminate early when
 * vertex 0 is won by player solvelocal. If solvelocal is -1 we are not allowed terminate early, if it is set to 2 we
 * are allowed to terminate when vertex 0 is found for either player 0 or 1.
 *
 * When solving globally solvelocal should be set to -1 initially.
 * Wehn solving locally solvelocal should be set to 2 initially.
 */
class zlnkPG {
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
     * Flag to enable metrics.
     *
     * Currently this does nothing but is present to keep the public variables inline with the public variables of
     * zlnkVPG.
     */
    static bool conf_metricoutput;
    /**
     * Measure time spent in the attractor set calculation
     */
    long attracting = 0;

    /**
     * Initialize solver
     * @param game Parity game
     */
    explicit zlnkPG(Game * game);
    /**
     * Initialize solver for the subgame containing vertices in bigV
     * @param game Parity game
     * @param bigV Vertices in the game
     */
    zlnkPG(Game * game, VertexSetZlnk * bigV);
    /**
     * Solve the parity game.
     * Assume the winning sets are empty.
     *
     * @param W0bigV Output parameter: winning set for player 0
     * @param W1bigV Output parameter: winning set for player 1
     */
    void solve(VertexSetZlnk * W0bigV, VertexSetZlnk * W1bigV);

protected:
    /**
     * Parity game
     */
    Game * game;
    /**
     * Indicate what subgame we are solving
     */
    VertexSetZlnk * bigV;

    /**
     * Attract vertices from bigA for player "player". Vertices attracted are removed from bigV, so attracting creates a
     * subgame simultaneously.
     * @param player Player for which to attract
     * @param bigA Vertices we are attracting
     */
    void attr(int player,VertexSetZlnk * bigA);
    /**
     * Attractor implementation using a queue to maintain which vertices were attracted and for we consider predecessors
     * Runs in O(e)
     */
    void attrQueue(int player,VertexSetZlnk * bigA);
    /**
     * Find the highest and lowest priority
     * @return tuple containing the lowest priority first and the highest second
     */
    tuple<int,int> getHighLowPrio();
    /**
     * Get all the vertices with a specific priority
     * @param bigA Output argument in which we place the vertices found
     * @param prio priority we are looking for
     */
    void getVCWithPrio(VertexSetZlnk *bigA, int prio);
    /**
     * Calculate the union of bigA and bigB and store the result in bigA
     * @param bigA Left hand set
     * @param bigB Right hand set
     */
    void unify(VertexSetZlnk * bigA, VertexSetZlnk * bigB);
    /**
     * Remove vertices in bigA from bigV
     * @param bigA Vertices to remove
     */
    void removeFromBigV(VertexSetZlnk * bigA);
    /**
     * Check if we are allowed to terminate early if we find vertex 0 for player "player"
     */
    bool inSolveLocal(int player);
};


#endif //VPGSOLVER_ZLNKPG_H
