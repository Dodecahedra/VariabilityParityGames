//
// Created by sjef on 11-7-19.
//

#ifndef VPGSOLVER_MBR_H
#define VPGSOLVER_MBR_H


#include "Game.h"
#include "conf.h"
#include "Measurements/bintree.h"
#include <vector>
#include <iostream>
using namespace std;

class MBR {
public:
    static vector<ConfSet> winningConf;
    static vector<VertexSetFPIte> winningVertices;
    Game * game;
    ConfSet * conf;
    VertexSetFPIte * P0;
    VertexSetFPIte * VP1;
    int feature;
    static bool fulloutput;

    bintree<vector<long>> *measured = nullptr;
    bool solvelocal = false;
    long confstring = 1;

    static bool metric_output;
    static string metric_dir;
    MBR(Game * game);
    MBR(Game * game, ConfSet * conf, VertexSetFPIte * P0, VertexSetFPIte * P1, int feature);

    void solve();
    void parition(ConfSet * org, ConfSet * part);
    void createSubGames(ConfSet * confP, vector<std::tuple<int,int>> * subgame_out, vector<std::tuple<int,int>> * subgame_in);
    void createPessimisticGames(vector<std::tuple<int,int>> * pessimistic_out0, vector<std::tuple<int,int>> *pessimistic_in0, vector<std::tuple<int,int>> *pessimistic_out1, vector<std::tuple<int,int>> *pessimistic_in1);
    void copyEdges(vector<std::tuple<int,int>> * edgeout, vector<std::tuple<int,int>> * edgein);
    void printMeasurements(ostream * output);
    int printNode(ostream * output,bintree<vector<long>> * node, int c);

    void removeDisabledEdge(vector<tuple<int, int>> *edge, vector<bool> *edgeenabled);
};


#endif //VPGSOLVER_MBR_H
