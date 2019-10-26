//
// Created by sjef on 26-10-19.
//

#ifndef VPGSOLVER_ZLNKPG_H
#define VPGSOLVER_ZLNKPG_H


#include "Game.h"

class zlnkPG {
public:
    static bool conf_metricoutput;
    static VertexSetZlnk emptyvertexset;
    long attracting = 0;
    Game * game;
    VertexSetZlnk * bigV;
    int solvelocal = -1;


    zlnkPG(Game * game);
    zlnkPG(Game * game, VertexSetZlnk * bigV);

    void attr(int player,VertexSetZlnk * bigA);
    void attrQueue(int player,VertexSetZlnk * bigA);
    void solve(VertexSetZlnk * W0bigV, VertexSetZlnk * W1bigV);
    tuple<int,int> getHighLowPrio();
    void getVCWithPrio(VertexSetZlnk *bigA, int prio);
    void unify(VertexSetZlnk * bigA, VertexSetZlnk * bigB);
    void removeFromBigV(VertexSetZlnk * bigA);
    void removeFromBigV(int i, ConfSet c);

    bool inSolveLocal(int player);
};


#endif //VPGSOLVER_ZLNKPG_H
