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
    vector<bdd> * vc;


    zlnk(Game * game);
    zlnk(Game * game, unordered_set<int> * bigV, vector<bdd> * vc);

    void attr(int player,unordered_set<int> * bigA,vector<bdd> * ac);
    void attrNaive(int player,unordered_set<int> * bigA,vector<bdd> * ac);
    void attrQueue(int player,unordered_set<int> * bigA,vector<bdd> * ac);
    void solve(unordered_set<int> * W0bigV, vector<bdd> * W0vc,unordered_set<int> * W1bigV, vector<bdd> * W1vc);
    tuple<int,int> getHighLowPrio();
    void getVCWithPrio(unordered_set<int> *bigA, vector<bdd> *ac, int prio);
    void unify(unordered_set<int> * bigA, vector<bdd> *ac, unordered_set<int> * bigB, vector<bdd> * bc);
    void removeFromBigV(unordered_set<int> * bigA, vector<bdd> *ac);
    void removeFromBigV(int i, bdd c);
    void test();
};


#endif //VPGSOLVER_ZLNK_H
