//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_ZLNK_H
#define VPGSOLVER_ZLNK_H


#include "Game.h"

class zlnk {
public:
    static bool conf_metricoutput;
    static VertexSetZlnk emptyvertexset;
    long attracting = 0;
    Game * game;
    VertexSetZlnk * bigV;
    vector<Subset> * vc;
    int solvelocal = -1;


    zlnk(Game * game);
    zlnk(Game * game, VertexSetZlnk * bigV, vector<Subset> * vc);

    void attr(int player,VertexSetZlnk * bigA,vector<Subset> * ac);
    void attrNaive(int player,VertexSetZlnk * bigA,vector<Subset> * ac);
    void attrQueue(int player,VertexSetZlnk * bigA,vector<Subset> * ac);
    void solve(VertexSetZlnk * W0bigV, vector<Subset> * W0vc,VertexSetZlnk * W1bigV, vector<Subset> * W1vc);
    tuple<int,int> getHighLowPrio();
    void getVCWithPrio(VertexSetZlnk *bigA, vector<Subset> *ac, int prio);
    void unify(VertexSetZlnk * bigA, vector<Subset> *ac, VertexSetZlnk * bigB, vector<Subset> * bc);
    void removeFromBigV(VertexSetZlnk * bigA, vector<Subset> *ac);
    void removeFromBigV(int i, Subset c);
    void test();

    void removeCF(Subset c, VertexSetZlnk * bigA, vector<Subset> *ac);
    bool inSolveLocal(int player);
};


#endif //VPGSOLVER_ZLNK_H
