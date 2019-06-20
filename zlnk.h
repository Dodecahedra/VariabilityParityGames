//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_ZLNK_H
#define VPGSOLVER_ZLNK_H


#include <unordered_set>
#include "Game.h"

class zlnk {
public:
    long attracting = 0;
    Game * game;
    unordered_set<int> * bigV;
    vector<Subset> * vc;


    zlnk(Game * game);
    zlnk(Game * game, unordered_set<int> * bigV, vector<Subset> * vc);

    void attr(int player,unordered_set<int> * bigA,vector<Subset> * ac);
    void attrNaive(int player,unordered_set<int> * bigA,vector<Subset> * ac);
    void attrQueue(int player,unordered_set<int> * bigA,vector<Subset> * ac);
    void solve(unordered_set<int> * W0bigV, vector<Subset> * W0vc,unordered_set<int> * W1bigV, vector<Subset> * W1vc);
    tuple<int,int> getHighLowPrio();
    void getVCWithPrio(unordered_set<int> *bigA, vector<Subset> *ac, int prio);
    void unify(unordered_set<int> * bigA, vector<Subset> *ac, unordered_set<int> * bigB, vector<Subset> * bc);
    void removeFromBigV(unordered_set<int> * bigA, vector<Subset> *ac);
    void removeFromBigV(int i, Subset c);
    void test();
};


#endif //VPGSOLVER_ZLNK_H
