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
    vector<int> successorsnotinbigA; // count the number of successors not attracted
    vector<bool> countinitialized; // indicate for every vertex if successorsnotinbigA is initialized
    successorsnotinbigA.resize(game->n_nodes);
    countinitialized.resize(game->n_nodes);
    queue<int> qq;
#ifdef VertexSetZlnkIsBitVector //@todo: implement some sort of iteration for VertexSetZlnkIsBitVector
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
            // Since everything in bigA is removed from bigA we conclude that vi is not in bigA

            // vertex vii has been placed into the attractor set, so vi has at least one successor that is attracted, namely vii
            // if vertex vi is owned by player "player" then we can attract vi
            bool attracted = true;
            if(game->owner[vi] != player){
                // if vi is not owned by player "player" then we need to check that all successors are in bigA
                if(!countinitialized[vi]){
                    countinitialized[vi] = true;
                    successorsnotinbigA[vi] = 0;
                    for(auto & j : game->out_edges[vi])
                        if((*bigV)[target(j)] || (*bigA)[target(j)])
                            successorsnotinbigA[vi]++;
                }
                // Maintain successorsnotinbigA value, only if there are no successors not in bigA can we attract vi
                if((--successorsnotinbigA[vi]) > 0){
                    attracted = false;
                }
            }
            if(!attracted)
                continue;
            // attracted vi, place it in the queue
            (*bigA)[vi] = true;
            (*bigV)[vi] = false;
            qq.push(vi);
        }
    }
}

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


    // initialize bigA which is used to attract
    auto * bigA = new VertexSetZlnk(game->n_nodes);
    // initialize subBigV which describes the subgame
    auto * subBigV = new VertexSetZlnk(game->n_nodes);
    // fill bigA with highest priorities
    getVCWithPrio(bigA, h);
    // initially, fill subgame with the entire game
    *subBigV = *bigV;
    //create subgame
    zlnkPG subgame(game, subBigV);
    // attract bigA, everything that is attracted is removed from subBigV, thus creating the correct subgame
    subgame.attr(player, bigA);
    cout << "\nDown1\n";

    // Everything won by player "1 - player" in the subgame is also won by "1 - player" in this game.
    // If we are allowed to terminate early when finding vertex 0 for "1 - player" then we are tell the subgame it is
    // also allowed to termiante early for "1 - player", otherwise it is never allowed to terminate early.
    // It is never allowed to terminate early when vertex 0 is found for player "player" since we are not sure if vertex
    // 0 is then also won by "player" in this game.
    if(inSolveLocal(1 - player))
        subgame.solvelocal = 1-player;
    else
        subgame.solvelocal = -1;

    subgame.solve(W0bigV, W1bigV);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    if(*WOpbigV == zlnkPG::emptyvertexset){
        // The entire subgame is won by player "player" and therefore every vertex in this game is won by "player"
        unify(WMebigV, bigA);
    } else {
        if(inSolveLocal(1-player) && (*WOpbigV)[0]){
            // Vertex 0 is won by player "1 - player" in the subgame and therefore also won by "1 - player" in this game
            // We don't have to go in the second recursion, we can return the winning sets as is.

            // Note: this weird notation (empty if statement) is used to keep the decision structure inline with the
            // pseudo code.
        } else {
            // clone content and wipe winningConf sets
            *bigA = *WOpbigV;
#ifdef VertexSetZlnkIsBitVector //@todo: implement some sort clear function VertexSetZlnkIsBitVector
            std::fill(W0bigV->begin(), W0bigV->end(), false);
            std::fill(W1bigV->begin(), W1bigV->end(), false);
#else
            W0bigV->clear();
            W1bigV->clear();
#endif
            // Attract bigA (which is now contains the vertices won by player "1- player" in the subgame) for player "1 - player"
            zlnkPG subgame2(game, bigV);
            subgame2.attr(1 - player, bigA);
            if(inSolveLocal(1-player) && (*bigA)[0]){
                // If vertex 0 is in the attracted set then player "1 - player" can force the play to a vertex that was
                // won by "1 - player" in the subgame, since these vertices are also won by "1 - player" in this game
                // we conclude vertex 0 is won by "1 - player" in this game and we can terminate.
                *WOpbigV = *bigA;
            } else {
                // Calculate final winning sets.
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
#ifdef VertexSetZlnkIsBitVector //@todo: implement some sort of iteration for VertexSetZlnkIsBitVector
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
    // use inverse priority assignment function
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