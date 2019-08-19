//
// Created by sjef on 5-6-19.
//

#include <queue>
#include <iostream>
#include "zlnk.h"
#include <map>
#include <chrono>


bool zlnk::conf_metricoutput = false;
VertexSetZlnk zlnk::emptyvertexset;
zlnk::zlnk(Game *game) {
    zlnk::emptyvertexset = VertexSetZlnk(game->n_nodes);
#ifdef SINGLEMODE
    vector<Subset> * vc = nullptr;
#else
    vector<Subset> * vc = new vector<Subset>(game->n_nodes);
#endif
    bigV = new VertexSetZlnk(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++){
        (*bigV)[i] = true;
#ifndef SINGLEMODE
        (*vc)[i] = game->bigC;
#endif
    }
    this->game = game;
    this->vc = vc;
}

zlnk::zlnk(Game * game, VertexSetZlnk * bigV, vector<Subset> * vc){
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

void zlnk::attr(int player, VertexSetZlnk *bigA, vector<Subset> *ac) {
    auto start = std::chrono::high_resolution_clock::now();



    removeFromBigV(bigA, ac);
    attrQueue(player, bigA, ac);
    auto end = std::chrono::high_resolution_clock::now();





    auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Attracting took " << elapsed.count() << "ns";
    attracting += elapsed.count();
}

#ifdef SINGLEMODE
void zlnk::attrQueue(int player, VertexSetZlnk *bigA, vector<Subset> *ac) {
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
#else
void zlnk::attrQueue(int player, VertexSetZlnk *bigA, vector<Subset> *ac) {
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
            int gi = guard_index(game->in_edges[vii][i]);
            if(!(*bigV)[vi]) // vertex not in the playing area anymore
                continue;
            Subset attracted;
            if(game->owner[vi] == player){
                attracted = (*vc)[vi];
                attracted &= (*ac)[vii];
                attracted &= game->edge_guards[gi];
            } else {
                attracted = (*vc)[vi];
                for(auto & j : game->out_edges[vi]){
                    int target = target(j);
                    Subset s = game->bigC;
                    Subset s2 = game->bigC;
                    s -= game->edge_guards[guard_index(j)];
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
            (*vc)[vi] -= attracted;
            if((*vc)[vi] == emptyset)
            {
                (*bigV)[vi] = false;
            }
            qq.push(vi);
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
#endif

#ifdef SINGLEMODE
/**
 * Assume empty winning sets (bigV empty and vc all 0's)
 * @param W0bigV
 * @param W0vc
 * @param W1bigV
 * @param W1vc
 */
void zlnk::solve(VertexSetZlnk *W0bigV, vector<Subset> *W0vc, VertexSetZlnk *W1bigV, vector<Subset> *W1vc) {
    if(*bigV == zlnk::emptyvertexset) return;

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
    getVCWithPrio(bigA, nullptr, h);

    *subBigV = *bigV;
    zlnk subgame(game, subBigV, nullptr);
    subgame.attr(player, bigA, nullptr);
    cout << "\nDown1\n";
    subgame.solve(W0bigV, W0vc, W1bigV, W1vc);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    if(*WOpbigV == zlnk::emptyvertexset){
        unify(WMebigV, nullptr, bigA, nullptr);
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

        zlnk subgame2(game, bigV, vc);
        subgame2.attr(1-player, bigA, nullptr);
        cout << "\nDown2\n";
        subgame2.solve(W0bigV, W0vc, W1bigV, W1vc);
        attracting += subgame.attracting;
        cout << "\nUp2\n";
        unify(WOpbigV, nullptr, bigA, nullptr);
    }
    delete bigA;
    delete subBigV;
}
#else
/**
 * Assume empty winningConf sets (bigV empty and vc all 0's)
 * @param W0bigV
 * @param W0vc
 * @param W1bigV
 * @param W1vc
 */
void zlnk::solve(VertexSetZlnk *W0bigV, vector<Subset> *W0vc, VertexSetZlnk *W1bigV, vector<Subset> *W1vc) {
    if(*bigV == zlnk::emptyvertexset) return;

    auto [h,l] = getHighLowPrio();
    int player = (h % 2);
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
    vector<Subset> *WMevc;
    VertexSetZlnk *WOpbigV;
    vector<Subset> *WOpvc;
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



    auto * bigA = new VertexSetZlnk(game->n_nodes);;
    vector<Subset> * ac = new vector<Subset>(game->n_nodes);

    auto * subBigV = new VertexSetZlnk(game->n_nodes);;
    vector<Subset> * subvc = new vector<Subset>(game->n_nodes);

    getVCWithPrio(bigA, ac, h);

    *subBigV = *bigV;
    *subvc = *vc;
    zlnk subgame(game, subBigV, subvc);
    subgame.attr(player, bigA, ac);
    cout << "\nDown1\n";
    subgame.solve(W0bigV, W0vc, W1bigV, W1vc);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    if(*WOpbigV == zlnk::emptyvertexset){
        unify(WMebigV, WMevc, bigA, ac);
    } else {
        // clone content and wipe winningConf sets
        *bigA = *WOpbigV;
        *ac = *WOpvc;
#ifdef VertexSetZlnkIsBitVector
        std::fill(W0bigV->begin(), W0bigV->end(), false);
        std::fill(W1bigV->begin(), W1bigV->end(), false);
#else
        W0bigV->clear();
        W1bigV->clear();
#endif
        fill(W0vc->begin(), W0vc->end(), emptyset);
        fill(W1vc->begin(), W1vc->end(), emptyset);

        zlnk subgame2(game, bigV, vc);
        subgame2.attr(1-player, bigA, ac);
        cout << "\nDown2\n";
        subgame2.solve(W0bigV, W0vc, W1bigV, W1vc);
        attracting += subgame.attracting;
        cout << "\nUp2\n";
        unify(WOpbigV, WOpvc, bigA, ac);
    }
    delete bigA;
    delete ac;
    delete subBigV;
    delete subvc;
}
#endif

tuple<int, int> zlnk::getHighLowPrio() {
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

void zlnk::getVCWithPrio(VertexSetZlnk *bigA, vector<Subset> *ac, int prio) {
    for (const auto& vi : game->priorityI[prio]) {
        if((*bigV)[vi]){
            (*bigA)[vi] = true;
#ifndef SINGLEMODE
            (*ac)[vi] = (*vc)[vi];
#endif
        }
    }
}

void zlnk::unify(VertexSetZlnk *bigA, vector<Subset> *ac, VertexSetZlnk *bigB, vector<Subset> *bc) {
#ifdef VertexSetZlnkIsBitVector
    for (int vi = 0;vi < game->n_nodes;vi++){
        if(!(*bigB)[vi])
            continue;
#else
    for (const auto& vi : *bigB) {
#endif
        (*bigA)[vi] = true;
#ifndef SINGLEMODE
        (*ac)[vi] |= (*bc)[vi];
#endif
    }
}

#ifdef SINGLEMODE
void zlnk::removeFromBigV(VertexSetZlnk * bigA, vector<Subset> *ac){
#ifdef VertexSetZlnkIsBitVector
    for (int vi = 0;vi < game->n_nodes;vi++){
        if(!(*bigA)[vi])
            continue;
#else
    for (const auto& vi : *bigA) {
#endif
        (*bigV)[vi] = false;
    }
}
void zlnk::removeFromBigV(int i, Subset c) {
    (*bigV)[i] = false;
}

#else
void zlnk::removeFromBigV(VertexSetZlnk * bigA, vector<Subset> *ac){
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

void zlnk::removeFromBigV(int i, Subset c) {
    (*vc)[i] -= c;
    if((*vc)[i] == emptyset)
    {
        (*bigV)[i] = false;
    }
}

#endif