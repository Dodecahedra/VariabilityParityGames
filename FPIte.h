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
    Game * game;
    VertexSet * P0;
    VertexSet * VP1;
    int d;
    VertexSet * W0;
    vector<bool> * edgeenabled;

    vector<VertexSet> ZZ;
    FPIte(Game * game);
    FPIte(Game * game, VertexSet * P0, VertexSet * VP1,vector<bool> * edgeenabled);

    void init(int i);
    void solve();
    void diamondbox(VertexSet * Z);
};


#endif //VPGSOLVER_FPITE_H
