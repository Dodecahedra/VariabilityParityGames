//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_ZLNK_H
#define VPGSOLVER_ZLNK_H


#include <unordered_set>
#include "Game.h"

class zlnk {
public:
    Game * game;
    unordered_set<int> * bigV;
    vector<BDD> * vc;


    zlnk(Game * game);
    zlnk(Game * game, unordered_set<int> * bigV, vector<BDD> * vc);

    void attr(int player, std::vector<tuple<int, BDD>> a, int ia);
};


#endif //VPGSOLVER_ZLNK_H
