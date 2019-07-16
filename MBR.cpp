//
// Created by sjef on 11-7-19.
//

#include "MBR.h"
#include "FPIte.h"
#include <iostream>
#include <algorithm>
vector<Subset> MBR::winningConf;
vector<VertexSet> MBR::winningVertices;

MBR::MBR(Game *game, Subset *conf, vector<bool> *edgeenabled, VertexSet *P0, VertexSet *VP1, int feature) {
    this->game = game;
    this->conf = conf;
    this->edgeenabled = edgeenabled;
    this->P0 = P0;
    this->VP1 = VP1;
    this->feature = feature;
}

MBR::MBR(Game *game) {
    this->game = game;
    this->conf = new Subset;
    *this->conf = game->bigC;

    this->P0 =  new VertexSet();
    this->P0->resize(game->n_nodes);
    this->VP1 =  new VertexSet();
    this->VP1->resize(game->n_nodes);
    fill(this->VP1->begin(), this->VP1->end(), true);
    this->edgeenabled = new vector<bool>(game->edge_guards.size());
    fill(this->edgeenabled->begin(), this->edgeenabled->end(), true);
    this->feature =  0;
}

void MBR::solve() {
    if(this->feature == game->bm_n_vars){
        int i = winningConf.size();
        winningConf.resize(i+1);
        winningVertices.resize(i+1);
        winningConf[i] = *conf;
        winningVertices[i].resize(game->n_nodes);
        FPIte fpite(game, P0, VP1, edgeenabled, &(winningVertices[i]));
        fpite.solve();
        return;
    }
#ifdef subsetbdd
    int sizeP0 = std::count_if(P0->begin(), P0->end(), [](bool b){return b;});
    int sizeP1 = game->n_nodes - std::count_if(VP1->begin(), VP1->end(), [](bool b){return b;});
    cout << "Size P0: " << sizeP0 << "/" << game->n_nodes << endl;
    cout << "Size P1: " << sizeP1 << "/" << game->n_nodes << endl;
#endif

    vector<bool> pessimisticedges0(edgeenabled->size());
    vector<bool> pessimisticedges1(edgeenabled->size());

    createPessimisticGames(&pessimisticedges0, &pessimisticedges1);

    VertexSet W0(game->n_nodes);
    auto * fpite0 = new FPIte(game, P0, VP1, &pessimisticedges0, &W0);
    fpite0->solve();
    *P0 = W0;
    auto * fpite1 = new FPIte(game, P0, VP1, &pessimisticedges1, &W0);
    fpite1->solve();
    *VP1 = W0;

    delete fpite0;
    delete fpite1;
    bool done = *P0 == *VP1;
    if(done){
        int i = winningConf.size();
        winningConf.resize(i+1);
        winningVertices.resize(i+1);
        winningConf[i] = *conf;
        winningVertices[i] = *P0;
        return;
    }

    auto * confa = conf;
    auto * confb = new Subset;
    parition(confa,confb);

    fill(pessimisticedges0.begin(), pessimisticedges0.end(), false);
    createSubGames(confa, &pessimisticedges0);
    createSubGames(confb, edgeenabled);

    auto * P0b = new VertexSet;
    auto * VP1b = new VertexSet;
    *P0b = *P0;
    *VP1b = *VP1;

//    fill(P0->begin(), P0->end(), false);
//    fill(P0b->begin(), P0b->end(), false);
//    fill(VP1->begin(), VP1->end(), true);
//    fill(VP1b->begin(), VP1b->end(), true);

    MBR ma(game, confa, &pessimisticedges0, P0, VP1, feature+1);
    MBR mb(game, confb, edgeenabled, P0b, VP1b, feature+1);
    ma.solve();
    mb.solve();

    delete P0b;
    delete VP1b;
}

void MBR::createPessimisticGames(vector<bool> *pessimisticedges0, vector<bool> *pessimisticedges1) {
    vector<bool> * pess;
    vector<bool> * opt;
    for(int i = 0;i<edgeenabled->size();i++){
        if(!(*edgeenabled)[i])
            continue;
        if(game->owner[game->edge_origins[i]] == 0){
            pess = pessimisticedges0;
            opt = pessimisticedges1;
        } else{
            pess = pessimisticedges1;
            opt = pessimisticedges0;
        }
        (*opt)[i] = true;
        Subset g = *conf;
        g -= game->edge_guards[i];
        (*pess)[i] = g == emptyset;
    }
}

void MBR::parition(Subset *org, Subset *part) {
    *part = *org;
    *org &= game->bm_vars[feature];
    *part -= game->bm_vars[feature];
}

void MBR::createSubGames(Subset *confP, vector<bool> *edgeenabledP) {
    for(int i = 0;i<edgeenabled->size();i++) {
        if (!(*edgeenabled)[i])
            continue;
        Subset g = *confP;
        g &= game->edge_guards[i];
        (*edgeenabledP)[i] = !(g == emptyset);
    }
}
