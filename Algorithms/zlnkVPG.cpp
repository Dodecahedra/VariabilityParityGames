//
// Created by sjef on 5-6-19.
//

#include "zlnkVPG.h"
#include <queue>
#include <iostream>
#include <map>
#include <chrono>


bool zlnkVPG::conf_metricoutput = false;
VertexSetZlnk zlnkVPG::emptyvertexset;

zlnkVPG::zlnkVPG(VPGame *game) {
    zlnkVPG::emptyvertexset = VertexSetZlnk(game->n_nodes);
    vector<ConfSet> * vc = new vector<ConfSet>(game->n_nodes);
    bigV = new VertexSetZlnk(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++){
        (*bigV)[i] = true;
        (*vc)[i] = game->bigC;
    }
    this->game = game;
    this->vc = vc;
}

zlnkVPG::zlnkVPG(VPGame * game, VertexSetZlnk * bigV, vector<ConfSet> * vc){
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

void zlnkVPG::attr(int player, VertexSetZlnk *bigA, vector<ConfSet> *ac) {
    auto start = std::chrono::high_resolution_clock::now();

    removeFromBigV(bigA, ac);
    attrQueue(player, bigA, ac);
    auto end = std::chrono::high_resolution_clock::now();

    auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Attracting took " << elapsed.count() << "ns";
    attracting += elapsed.count();
}


void zlnkVPG::attrQueue(int player, VertexSetZlnk *bigA, vector<ConfSet> *ac) {
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

        for(int i = 0;i<game->in_edges[vii].size();i++){
            int vi = target(game->in_edges[vii][i]);
            int gi = edge_index(game->in_edges[vii][i]);
            if(!(*bigV)[vi]) // vertex not in the playing area anymore
                continue;
            // try to attract as many configurations as possible for vertex vi
            ConfSet attracted;
            // This follows the attractor definition precisely, see preudo code and definitions in the report
            if(game->owner[vi] == player){
                attracted = (*vc)[vi];
                attracted &= (*ac)[vii];
                attracted &= game->edge_guards[gi];
            } else {
                attracted = (*vc)[vi];
                for(auto & j : game->out_edges[vi]){
                    int target = target(j);
                    ConfSet s = game->bigC;
                    ConfSet s2 = game->bigC;
                    s -= game->edge_guards[edge_index(j)];
                    s2 -= (*vc)[target];
                    s |= s2;
                    attracted &= s;
                }
            }
            if(attracted == emptyset)
                continue;
            if(!(*bigA)[vi]){
                // add vertex to attracted set
                (*bigA)[vi] = true;
                (*ac)[vi] = attracted;
            } else {
                (*ac)[vi] |= attracted;
            }
            // remove attracted confs from the game
            (*vc)[vi] -= attracted;
            if((*vc)[vi] == emptyset)
            {
                (*bigV)[vi] = false;
            }
            // if we attracted anything we need to reevaluate the predecessors of vi
            qq.push(vi);

            // Count how many configurations we were able to attract
            if(conf_metricoutput)
#ifdef subsetbdd
                cout << "Attracted " << bdd_satcount(attracted) << " configurations.\n";
#elif subsetexplicit
                cout << "Attracted " << attracted.count() << " configurations.\n";
#else
                ;
#endif
        }
    }
}

void zlnkVPG::solve(VertexSetZlnk *W0bigV, vector<ConfSet> *W0vc, VertexSetZlnk *W1bigV, vector<ConfSet> *W1vc) {
    if(*bigV == zlnkVPG::emptyvertexset) return;

    auto [h,l] = getHighLowPrio();
    int player = (l % 2);
    if(h == l){
        if(player == 0){
            *W0bigV = *bigV;
            *W0vc = *vc;
        } else {
            *W1bigV = *bigV;
            *W1vc = *vc;
        }
        return;
    }

    VertexSetZlnk *WMebigV;
    vector<ConfSet> *WMevc;
    VertexSetZlnk *WOpbigV;
    vector<ConfSet> *WOpvc;
    if(player == 0){
        WMebigV = W0bigV;
        WMevc = W0vc;
        WOpbigV = W1bigV;
        WOpvc = W1vc;
    } else {
        WMebigV = W1bigV;
        WMevc = W1vc;
        WOpbigV = W0bigV;
        WOpvc = W0vc;
    }



    // initialize (bigA,ac) which is used to attract
    auto * bigA = new VertexSetZlnk(game->n_nodes);;
    vector<ConfSet> * ac = new vector<ConfSet>(game->n_nodes);

    // initialize (subBigV,subvc) which describes the subgame
    auto * subBigV = new VertexSetZlnk(game->n_nodes);;
    vector<ConfSet> * subvc = new vector<ConfSet>(game->n_nodes);
    //fill (bigA,ac) with the lowest priorities
    getVCWithPrio(bigA, ac, l);
    // initially, fill subgame with the entire game
    *subBigV = *bigV;
    *subvc = *vc;
    //create subgame
    zlnkVPG subgame(game, subBigV, subvc);
    //attract (bigA, ac), everything that is attracted is removed from (subBigV, subvc) thus creating the correct subgame
    subgame.attr(player, bigA, ac);
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

    subgame.solve(W0bigV, W0vc, W1bigV, W1vc);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    if(*WOpbigV == zlnkVPG::emptyvertexset){
        // The entire subgame is won by player "player" for all configurations
        // therefore every vertex in this game is won by "player" for all configurations
        unify(WMebigV, WMevc, bigA, ac);
    } else {
        ConfSet localconfs, localconfs2;
        if(inSolveLocal(1 - player) && (*WOpbigV)[0])
            localconfs = (*WOpvc)[0]; // Vertex 0 is found for these configurations, we are done for these confs
        if(!(localconfs == emptyset) && !removeCSet(WOpbigV, WOpvc, localconfs)){ // remove all confs in localconfs from the game
            // WOp is empty
            // The entire subgame is won by player "player" for all configuration except those in localconfs
            // therefore every vertex in this game is won by "player" for all configurations not in localconfs

            // create incomplete winning sets for confs in localconfs
            (*WOpbigV)[0] = true;
            (*WOpvc)[0] |= localconfs;
            // other confs are winning for player "player"
            unify(WMebigV,WMevc,bigA,ac);
        } else {
            // clone content and wipe winningConf sets
            *bigA = *WOpbigV;
            *ac = *WOpvc;
    #ifdef VertexSetZlnkIsBitVector //@todo: implement some sort clear function VertexSetZlnkIsBitVector
            std::fill(W0bigV->begin(), W0bigV->end(), false);
            std::fill(W1bigV->begin(), W1bigV->end(), false);
    #else
            W0bigV->clear();
            W1bigV->clear();
    #endif
            fill(W0vc->begin(), W0vc->end(), emptyset);
            fill(W1vc->begin(), W1vc->end(), emptyset);

            // Attract (bigA,ac) (which is now contains the vertices won by player "1- player" in the subgame) for
            // player "1 - player"
            zlnkVPG subgame2(game, bigV, vc);
            subgame2.attr(1-player, bigA, ac);
            if(inSolveLocal(1 - player)){
                localconfs2 = (*ac)[0]; // Vertex 0 is found for these configurations, we are done for these confs
                localconfs |= localconfs2;
            }
            // remove the newly found confs from the game and only continue of something is left
            if(!inSolveLocal(1 - player) || removeCSet(bigV, vc, localconfs2)){
                // Go into the second recursion for the remaining confs
                cout << "\nDown2\n";
                subgame2.solvelocal = solvelocal;
                subgame2.solve(W0bigV, W0vc, W1bigV, W1vc);
                attracting += subgame.attracting;
                cout << "\nUp2\n";
                unify(WOpbigV, WOpvc, bigA, ac);
            }
            if(!(localconfs == emptyset)) {
                // Add incomplete winning sets for the localconfs found
                (*WOpbigV)[0] = true;
                (*WOpvc)[0] |= localconfs;
            }
        }
    }
    delete bigA;
    delete ac;
    delete subBigV;
    delete subvc;
}

tuple<int, int> zlnkVPG::getHighLowPrio() {
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

void zlnkVPG::getVCWithPrio(VertexSetZlnk *bigA, vector<ConfSet> *ac, int prio) {
    // use inverse priority assignment function
    for (const auto& vi : game->priorityI[prio]) {
        if((*bigV)[vi]){
            (*bigA)[vi] = true;
            (*ac)[vi] = (*vc)[vi];
        }
    }
}

void zlnkVPG::unify(VertexSetZlnk *bigA, vector<ConfSet> *ac, VertexSetZlnk *bigB, vector<ConfSet> *bc) {
#ifdef VertexSetZlnkIsBitVector
    for (int vi = 0;vi < game->n_nodes;vi++){
        if(!(*bigB)[vi])
            continue;
#else
    for (const auto& vi : *bigB) {
#endif
        (*bigA)[vi] = true;
        (*ac)[vi] |= (*bc)[vi];
    }
}

void zlnkVPG::removeFromBigV(VertexSetZlnk * bigA, vector<ConfSet> *ac){
#ifdef VertexSetZlnkIsBitVector
    for (int vi = 0;vi < game->n_nodes;vi++){
        if(!(*bigA)[vi])
            continue;
#else
        for (const auto& vi : *bigA) {
#endif
        removeFromBigV(vi, (*ac)[vi]);
    }
}

void zlnkVPG::removeFromBigV(int i, ConfSet c) {
    (*vc)[i] -= c;
    if((*vc)[i] == emptyset)
    {
        (*bigV)[i] = false;
    }
}


bool zlnkVPG::inSolveLocal(int player) {
    if(solvelocal == 2)
        return true;
    return player == solvelocal;
}

bool zlnkVPG::removeCSet(VertexSetZlnk *bigA, vector<ConfSet> *ac, ConfSet C) {
    if(C == emptyset)
        return true;
    bool somethingleft = false;
    for (int vi = 0;vi < game->n_nodes;vi++) {
        if (!(*bigA)[vi])
            continue;
        (*ac)[vi] -= C;
        if((*ac)[vi] == emptyset)
            (*bigA)[vi] = false;
        else
            somethingleft = true;
    }
    return somethingleft;
}
