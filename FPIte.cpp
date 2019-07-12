//
// Created by sjef on 11-7-19.
//

#include "FPIte.h"

#define targetIsIn(t) ZZ[game->priority[t]][t]
FPIte::FPIte(Game *game, VertexSet *P0, VertexSet *VP1, vector<bool> *edgeenabled, VertexSet * W0) {
    this->game = game;
    this->P0 = P0;
    this->VP1 = VP1;
    this->edgeenabled = edgeenabled;
    this->W0 = W0;
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
    this->W0 = new VertexSet();
    this->W0->resize(game->n_nodes);
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


    for(int i = 0;i<d;i++){
        init(i);
    }

    vector<VertexSet> ZZa;
    ZZa.resize(d);
    for(int i = 0;i<d;i++)
        ZZa[i].resize(game->n_nodes);
    int i;
    do {
        copyWithPrio(&ZZa[0],&ZZ[0],0);
        diamondbox(W0);
        copyWithPrio(&ZZ[0],W0,0);
        i = 0;
        while(compareWithPrio(&ZZ[i], &ZZa[i],i) && i < d -1){
            i++;
            copyWithPrio(&ZZa[i],&ZZ[i],i);
            copyWithPrio(&ZZ[i], W0,i);
            init(i-1);
        }
    } while (!(i == d-1 && compareWithPrio(&ZZ[d-1],&ZZa[d-1],d-1)));
}

FPIte::~FPIte() {
}

void FPIte::copyWithPrio(VertexSet *Z, VertexSet *ZP, int p) {
    for(auto& v : game->priorityI[p])
        (*Z)[v] = (*ZP)[v];
}

bool FPIte::copyAndCompareWithPrio(VertexSet *Z, VertexSet *ZP, int p) {
    bool equal = true;
    for(auto& v : game->priorityI[p]) {
        if((*Z)[v] != (*ZP)[v]) {
            equal = false;
            (*Z)[v] = (*ZP)[v];
        }
    }
    return equal;
}

bool FPIte::compareWithPrio(VertexSet *Z, VertexSet *ZP, int p) {
    bool equal = true;
    for(auto& v : game->priorityI[p])
        if((*Z)[v] != (*ZP)[v])
            equal = false;
    return equal;
}
