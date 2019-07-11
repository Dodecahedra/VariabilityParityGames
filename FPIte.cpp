//
// Created by sjef on 11-7-19.
//

#include "FPIte.h"

#define targetIsIn(t) ZZ[game->priority[t]][t]
FPIte::FPIte(Game *game, VertexSet *P0, VertexSet *VP1, vector<bool> *edgeenabled) {
    this->game = game;
    this->P0 = P0;
    this->VP1 = VP1;
    this->edgeenabled = edgeenabled;
}

FPIte::FPIte(Game *game) {
    this->game = game;
    this->P0 =  new VertexSet();
    this->P0->resize(game->n_nodes);
    this->VP1 =  new VertexSet();
    this->VP1->resize(game->n_nodes);
    fill(this->VP1->begin(), this->VP1->end(), true);
    this->edgeenabled = new vector<bool>(game->edge_guards.size());
    fill(this->edgeenabled->begin(), this->edgeenabled->end(), true);
}

void FPIte::init(int i) {
    if(i % 2 == 0)
        ZZ[i] = *VP1;
    else
        ZZ[i] = *P0;
}

void FPIte::diamondbox(VertexSet *Z) {
    for(int i = 0;i<game->n_nodes;i++){
        bool in;
        if(game->owner[i] == 0){
            in = false;
            for(int j = 0;j<game->out_edges[i].size() && !in;j++){
                auto edge = game->out_edges[i][j];
                if(!(*edgeenabled)[guard_index(edge)])
                    continue;
                int t = target(edge);
                if(targetIsIn(t))
                    in = true;
            }
        } else {
            in = true;
            for(int j = 0;j<game->out_edges[i].size() && in;j++) {
                auto edge = game->out_edges[i][j];
                if(!(*edgeenabled)[guard_index(edge)])
                    continue;
                int t = target(edge);
                if(!targetIsIn(t))
                    in = false;
            }
        }
        (*Z)[i] = in;
    }
}

void FPIte::solve() {
    d = game->priorityI.size();
    ZZ.resize(d);
    W0 = &(ZZ[d-1]);


    for(int i = 0;i<d;i++){
        init(i);
    }

    vector<VertexSet> ZZa;
    ZZa.resize(d);
    int i;
    do {
        ZZa[0] = ZZ[0];
        diamondbox(&ZZ[0]);
        i = 0;
        while(ZZ[i] == ZZa[i] && i < d -1){
            i++;
            ZZa[i] = ZZ[i];
            ZZ[i] = ZZ[i-1];
            init(i-1);
        }
    } while (!(i == d-1 && ZZ[d-1] == ZZa[d-1]));
}

FPIte::~FPIte() {
}
