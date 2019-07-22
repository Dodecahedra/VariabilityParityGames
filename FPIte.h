//
// Created by sjef on 11-7-19.
//

#ifndef VPGSOLVER_FPITE_H
#define VPGSOLVER_FPITE_H

#include "conf.h"
#include "Game.h"
#include <vector>
using namespace std;

class FPIte {
public:
    bool solvelocal = false;
    long verticesconsidered = 0;
    long dbs_executed = 0;
    Game * game;
    VertexSet * P0;
    VertexSet * VP1;
    int d;
    VertexSet * W0;
    vector<bool> * edgeenabled;
    vector<int> edgecount;

    VertexSet ZZ;
    VertexSet targetWasIn;
    FPIte(Game * game);
    FPIte(Game * game, VertexSet * P0, VertexSet * VP1,vector<bool> * edgeenabled, VertexSet * W0);
    ~FPIte();

    void init(int is, int ie);
    void solve();
    void diamondbox(VertexSet * Z);
    void diamondbox(VertexSet * Z, int maxprio);
    void copyWithPrio(VertexSet * Z, VertexSet * ZP, int p);
    void copyWithPrio(VertexSet * Z, VertexSet * ZP, int sp, int ep);
    bool copyAndCompareWithPrio(VertexSet * Z, VertexSet * ZP, int p);
    bool compareWithPrio(VertexSet * Z, VertexSet * ZP, int p);

    void setP0(char * P0string);
    void P0IsFull();
    void setP1(char * P1string);

    void unassist(int n);
};


#endif //VPGSOLVER_FPITE_H
