//
// Created by sjef on 11-7-19.
//

#ifndef VPGSOLVER_MBR_H
#define VPGSOLVER_MBR_H


#include "Game.h"
#include "conf.h"
#include "Algorithms/Datastructures/bintree.h"
#include <vector>
#include <iostream>
using namespace std;

/**
 * Implementation of the incremental pre-solve algorithm (working name is MBR)
 *
 * FPIte is used to solve (pessimistic) parity games
 *
 * @todo: experiment with different parity game algorithms to solve (pessimistic) parity games
 * @todo: make winning sets not static
 */
class MBR {
public:
    /**
     * Winning set for player 0, static because it is filled throughout the different recursions
     */
    static vector<ConfSet> winningConf;
    static vector<VertexSetFPIte> winningVertices;
    /**
     * Indicate if we want the fulloutput in the winning sets or if vertex 0 suffices
     */
    static bool fulloutput;
    /**
     * Flag to indicate if we are measuring stuff, this maintains a binary tree with the size of P0, VP1 and the time it
     * took to solve the pessimistic games
     */
    static bool metric_output;
    /**
     * Save all the pessimistic game created in this directory
     */
    static string metric_dir;
    /**
     * Indicate if we are allowed to terminate early when vertex 0 is found
     */
    bool solvelocal = false;

    /**
     * Initialize the solver
     * @param game Variability Parity game
     */
    explicit MBR(Game * game);
    /**
     * Initialize the solver. Solve only for configurations in conf and use P0 and P1. Create the next subgames by
     * splitting the set of configurations in two sets; one including feature "feature" and one excluding it
     * @param game VPG
     * @param conf configurations for which we are solving
     * @param P0 Vertices for which the VPG is won by player 0 for all configurations in conf
     * @param VP1 Vertices, such that every vertex not in VP1 is won by player 1 for all configurations in conf
     * @param feature Feature number on which we will split the set of configurations next
     */
    MBR(Game * game, ConfSet * conf, VertexSetFPIte * P0, VertexSetFPIte * VP1, int feature);

    /**
     * Solve the game, filling winningConf and winningVertices
     */
    void solve();

    /**
     * Print the binary tree measured in .dot format
     */
    void printMeasurements(ostream * output);
protected:
    /**
     * VPG
     */
    Game * game;
    /**
     * Describes the subgame
     */
    ConfSet * conf;
    /**
     * Vertices for which the VPG is won by player 0 for all configurations in conf
     */
    VertexSetFPIte * P0;
    /**
     * Vertices, such that every vertex not in VP1 is won by player 1 for all configurations in conf
     */
    VertexSetFPIte * VP1;
    /**
     * Feature number on which we will split the set of configurations next
     */
    int feature;

    /**
     * Binary tree measuring the size of P0, P1 and the solve times of the pessimistic parity games
     */
    bintree<vector<long>> *measured = nullptr;
    /**
     * String describing the set of configurations we are considering
     */
    long confstring = 1;

    /**
     * Parition the set of configuration, org is modified to include feature "feature" and part is modified to exclude
     * feature "feature"
     * @param org Current set of configurations
     * @param part New set of configurations
     */
    void parition(ConfSet * org, ConfSet * part);
    /**
     * Create new edge relation only containing configurations in confP, that is edges that do not contain any configuration
     * in confP are removed from the edge relation.
     * @param confP Set of configurations remaining
     * @param subgame_out Outgoing edge relation
     * @param subgame_in Incoming edge relation
     */
    void createSubGames(ConfSet * confP, vector<std::tuple<int,int>> * subgame_out, vector<std::tuple<int,int>> * subgame_in);
    /**
     * Create edge relations for the pessimistic parity games for player 0 and 1.
     * @param pessimistic_out0 Outgoing edge relation for the pessimistic parity game for player 0
     * @param pessimistic_in0 Incoming edge relation for the pessimistic parity game for player 0
     * @param pessimistic_out1 Outgoing edge relation for the pessimistic parity game for player 1
     * @param pessimistic_in1 Incoming edge relation for the pessimistic parity game for player 1
     */
    void createPessimisticGames(vector<std::tuple<int,int>> * pessimistic_out0, vector<std::tuple<int,int>> *pessimistic_in0, vector<std::tuple<int,int>> *pessimistic_out1, vector<std::tuple<int,int>> *pessimistic_in1);
    /**
     * Copy edge relation currently in game into edgeout and edgein
     * @param edgeout Outgoing edge relation in which the current outgoing edge relation is copied
     * @param edgein Incoming edge relation in which the current incoming edge relation is copied
     */
    void copyEdges(vector<std::tuple<int,int>> * edgeout, vector<std::tuple<int,int>> * edgein);
    /**
     * Remove edges from edge relation edge that are not enabled in edgeenabled
     * @param edge Edge relation from which to remove
     * @param edgeenabled Bool vector indicating if the edges should stay in or be removed
     */
    void removeDisabledEdge(vector<tuple<int, int>> *edge, vector<bool> *edgeenabled);


    /**
     * Print node in the binary tree to ostream
     */
    int printNode(ostream * output,bintree<vector<long>> * node, int c);
};


#endif //VPGSOLVER_MBR_H
