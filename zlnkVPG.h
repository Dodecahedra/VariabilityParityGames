//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_ZLNKVPG_H
#define VPGSOLVER_ZLNKVPG_H


#include "Game.h"

class zlnkVPG {
public:
    static bool conf_metricoutput;
    static VertexSetZlnk emptyvertexset;
    long attracting = 0;
    Game * game;
    VertexSetZlnk * bigV;
    vector<ConfSet> * vc;
    int solvelocal = -1;


    zlnkVPG(Game * game);
    zlnkVPG(Game * game, VertexSetZlnk * bigV, vector<ConfSet> * vc);

    void attr(int player,VertexSetZlnk * bigA,vector<ConfSet> * ac);
    void attrNaive(int player,VertexSetZlnk * bigA,vector<ConfSet> * ac);
    void attrQueue(int player,VertexSetZlnk * bigA,vector<ConfSet> * ac);
    void solve(VertexSetZlnk * W0bigV, vector<ConfSet> * W0vc,VertexSetZlnk * W1bigV, vector<ConfSet> * W1vc);
    tuple<int,int> getHighLowPrio();
    void getVCWithPrio(VertexSetZlnk *bigA, vector<ConfSet> *ac, int prio);
    void unify(VertexSetZlnk * bigA, vector<ConfSet> *ac, VertexSetZlnk * bigB, vector<ConfSet> * bc);
    void removeFromBigV(VertexSetZlnk * bigA, vector<ConfSet> *ac);
    void removeFromBigV(int i, ConfSet c);
    void test();

    void removeCF(ConfSet c, VertexSetZlnk * bigA, vector<ConfSet> *ac);
    bool inSolveLocal(int player);
    bool removeCSet(VertexSetZlnk *bigA, vector<ConfSet> *ac, ConfSet C);
};


#endif //VPGSOLVER_ZLNKVPG_H
