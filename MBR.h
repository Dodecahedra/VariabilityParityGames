//
// Created by sjef on 11-7-19.
//

#ifndef VPGSOLVER_MBR_H
#define VPGSOLVER_MBR_H


#include "Game.h"
#include "conf.h"
#include <vector>
using namespace std;

class MBR {
public:
    static vector<tuple<Subset, VertexSet>> winning;
    Game * game;
    Subset * conf;
    vector<bool> * edgeenabled;
    VertexSet * P0;
    VertexSet * P1;

    MBR(Game * game);
    MBR(Game * game, Subset * conf, vector<bool> * edgeenabled, VertexSet * P0, VertexSet * P1);

    void solve();
    void parition(Subset * org, Subset * part);
};


#endif //VPGSOLVER_MBR_H
