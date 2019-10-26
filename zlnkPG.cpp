//
// Created by sjef on 26-10-19.
//

#include "zlnkPG.h"
#include <queue>
#include <iostream>
#include <map>
#include <chrono>


bool zlnkPG::conf_metricoutput = false;
VertexSetZlnk zlnkPG::emptyvertexset;
zlnkPG::zlnkPG(Game *game) {
    zlnkPG::emptyvertexset = VertexSetZlnk(game->n_nodes);
    bigV = new VertexSetZlnk(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++){
        (*bigV)[i] = true;
    }
    this->game = game;
}

zlnkPG::zlnkPG(Game * game, VertexSetZlnk * bigV){
    this->game = game;
    this->bigV = bigV;
}

void zlnkPG::attr(int player, VertexSetZlnk *bigA) {
    auto start = std::chrono::high_resolution_clock::now();

    removeFromBigV(bigA);
    attrQueue(player, bigA);
    auto end = std::chrono::high_resolution_clock::now();



    auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Attracting took " << elapsed.count() << "ns";
    attracting += elapsed.count();
}

void zlnkPG::attrQueue(int player, VertexSetZlnk *bigA) {
    vector<bool> countinitialized;
    vector<int> countincoming;
    countincoming.resize(game->n_nodes);
    countinitialized.resize(game->n_nodes);
    queue<int> qq;
#ifdef VertexSetZlnkIsBitVector
    for(int vi = 0;vi<game->n_nodes;vi++){
        if(!(*bigA)[vi])
            continue;
#else
        for (const auto& vi : *bigA) {
#endif
        qq.push(vi);
    }
    while(!qq.empty())
    {
        int vii = qq.front();
        qq.pop();
        for(auto & i : game->in_edges[vii]){
            int vi = target(i);
            if(!(*bigV)[vi]) // vertex not in the playing area anymore
                continue;
            bool attracted = true;
            if(game->owner[vi] != player){
                if(!countinitialized[vi]){
                    countinitialized[vi] = true;
                    countincoming[vi] = 0;
                    for(auto & j : game->out_edges[vi])
                        if((*bigV)[target(j)] || (*bigA)[target(j)])
                            countincoming[vi]++;
                }
                if((--countincoming[vi]) > 0){
                    attracted = false;
                }
            }
            if(!attracted)
                continue;
            (*bigA)[vi] = true;
            (*bigV)[vi] = false;
            qq.push(vi);
        }
    }
}

/**
 * Assume empty winning sets (bigV empty and vc all 0's)
 * @param W0bigV
 * @param W0vc
 * @param W1bigV
 * @param W1vc
 */
void zlnkPG::solve(VertexSetZlnk *W0bigV, VertexSetZlnk *W1bigV) {
    cout << "Delta is " << solvelocal << endl;
    if(*bigV == zlnkPG::emptyvertexset) return;

    auto [h,l] = getHighLowPrio();
    int player = (h % 2);
    if(h == l){
        if(player == 0){
            *W0bigV = *bigV;
        } else {
            *W1bigV = *bigV;
        }
        return;
    }

    VertexSetZlnk *WMebigV;
    VertexSetZlnk *WOpbigV;
    if(player == 0){
        WMebigV = W0bigV;
        WOpbigV = W1bigV;
    } else {
        WMebigV = W1bigV;
        WOpbigV = W0bigV;
    }



    auto * bigA = new VertexSetZlnk(game->n_nodes);

    auto * subBigV = new VertexSetZlnk(game->n_nodes);
    getVCWithPrio(bigA, h);

    *subBigV = *bigV;
    zlnkPG subgame(game, subBigV);
    subgame.attr(player, bigA);
    cout << "\nDown1\n";
    if(inSolveLocal(1 - player))
        subgame.solvelocal = 1-player;
    else
        subgame.solvelocal = -1;

    subgame.solve(W0bigV, W1bigV);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    if(*WOpbigV == zlnkPG::emptyvertexset){
        unify(WMebigV, bigA);
    } else {
        if(inSolveLocal(1-player) && (*WOpbigV)[0]){
        } else {
            // clone content and wipe winningConf sets
            *bigA = *WOpbigV;
#ifdef VertexSetZlnkIsBitVector
            std::fill(W0bigV->begin(), W0bigV->end(), false);
            std::fill(W1bigV->begin(), W1bigV->end(), false);
#else
            W0bigV->clear();
            W1bigV->clear();
#endif

            zlnkPG subgame2(game, bigV);
            subgame2.attr(1 - player, bigA);
            if(inSolveLocal(1-player) && (*bigA)[0]){
                *WOpbigV = *bigA;
            } else {
                cout << "\nDown2\n";
                subgame2.solvelocal = solvelocal;
                subgame2.solve(W0bigV, W1bigV);
                attracting += subgame.attracting;
                cout << "\nUp2\n";
                unify(WOpbigV, bigA);
            }
        }
    }
    delete bigA;
    delete subBigV;
}


tuple<int, int> zlnkPG::getHighLowPrio() {
    int highest = 0;
    int lowest = INT32_MAX;
#ifdef VertexSetZlnkIsBitVector
    for(int vi = 0;vi<game->n_nodes;vi++){
        if(!(*bigV)[vi])
            continue;
#else
        for (const auto& vi : *bigV) {
#endif
        int prio = game->priority[vi];
        if(prio < lowest) lowest = prio;
        if(prio > highest) highest = prio;
    }
    return make_tuple(highest, lowest);
}

void zlnkPG::getVCWithPrio(VertexSetZlnk *bigA, int prio) {
    for (const auto& vi : game->priorityI[prio]) {
        if((*bigV)[vi]){
            (*bigA)[vi] = true;
        }
    }
}

void zlnkPG::unify(VertexSetZlnk *bigA, VertexSetZlnk *bigB) {
#ifdef VertexSetZlnkIsBitVector
    (*bigA) |= (*bigB);
#else
    for (const auto& vi : *bigB) {
        (*bigA)[vi] = true;
    }
#endif
}


void zlnkPG::removeFromBigV(VertexSetZlnk * bigA){
#ifdef VertexSetZlnkIsBitVector
    (*bigV) -= (*bigA);
//    for (int vi = 0;vi < game->n_nodes;vi++){
//        if(!(*bigA)[vi])
//            continue;
#else
    for (const auto& vi : *bigA) {
        (*bigV)[vi] = false;
    }
#endif
}

bool zlnkPG::inSolveLocal(int player) {
    if(solvelocal == 2)
        return true;
    return player == solvelocal;
}