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
    static vector<Subset> winningConf;
    static vector<VertexSet> winningVertices;
    Game * game;
    Subset * conf;
    vector<bool> * edgeenabled;
    VertexSet * P0;
    VertexSet * VP1;
    int feature;

    MBR(Game * game);
    MBR(Game * game, Subset * conf, vector<bool> * edgeenabled, VertexSet * P0, VertexSet * P1, int feature);

    void solve();
    void parition(Subset * org, Subset * part);
    void createSubGames(Subset * confP, vector<bool> * edgeenabledP);
    void createPessimisticGames(vector<bool> * pessimisticedges0, vector<bool> * pessimisticedges1);
};


#endif //VPGSOLVER_MBR_H
