//
// Created by sjef on 11-7-19.
//

#include "MBR.h"

vector<tuple<Subset, VertexSet>> MBR::winning;

MBR::MBR(Game *game, Subset *conf, vector<bool> *edgeenabled, VertexSet *P0, VertexSet *P1) {
    this->game = game;
    this->conf = conf;
    this->edgeenabled = edgeenabled;
    this->P0 = P0;
    this->P1 = P1;
}

MBR::MBR(Game *game) {
    this->game = game;
    this->conf = new Subset;
    *this->conf = game->bigC;

}
