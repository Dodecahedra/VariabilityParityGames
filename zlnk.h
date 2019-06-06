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

    void attr(int player,unordered_set<int> * bigA,vector<BDD> * ac);
    void attrNaive(int player,unordered_set<int> * bigA,vector<BDD> * ac);
    void attrQueue(int player,unordered_set<int> * bigA,vector<BDD> * ac);
    void solve(unordered_set<int> * W0bigV, vector<BDD> * W0vc,unordered_set<int> * W1bigV, vector<BDD> * W1vc);
    tuple<int,int> getHighLowPrio();
    void getVCWithPrio(unordered_set<int> *bigA, vector<BDD> *ac, int prio);
    void unify(unordered_set<int> * bigA, vector<BDD> *ac, unordered_set<int> * bigB, vector<BDD> * bc);
    void removeFromBigV(unordered_set<int> * bigA, vector<BDD> *ac);
    void removeFromBigV(int i, BDD c);
    void test();
};


#endif //VPGSOLVER_ZLNK_H
