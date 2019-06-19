//
// Created by sjef on 5-6-19.
//

#include <queue>
#include "zlnk.h"
#include <chrono>

zlnk::zlnk(Game *game) {
#ifdef SINGLEMODE
    vector<Subset> * vc = nullptr;
#else
    vector<Subset> * vc = new vector<Subset>(game->n_nodes);
#endif
    unordered_set<int> * bigV = new unordered_set<int>(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++){
        bigV->insert(i);
#ifndef SINGLEMODE
        (*vc)[i] = game->bigC;
#endif
    }
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

zlnk::zlnk(Game * game, unordered_set<int> * bigV, vector<Subset> * vc){
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

void zlnk::attr(int player, unordered_set<int> *bigA, vector<Subset> *ac) {
    auto start = std::chrono::high_resolution_clock::now();



    removeFromBigV(bigA, ac);
    attrQueue(player, bigA, ac);
    auto end = std::chrono::system_clock::now();





    auto elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Attracting took " << elapsed.count() << "ms.\n";
    attracting += elapsed.count();
}

#ifdef SINGLEMODE
void zlnk::attrQueue(int player, unordered_set<int> *bigA, vector<Subset> *ac) {
//    cout << "Content of bigV: ";
//    for(const auto& vi : *bigV){
//        cout << vi << ",";
//    }
//    cout << "\n";
    cout << "Attr start, size: " << bigA->size() << ", player: " << player << "\n";
    queue<int> qq;
    for (const auto& vi : *bigA) {
        qq.push(vi);
//        cout << vi << ",";
    }
    cout << "\nAttracted: ";
    while(!qq.empty())
    {
        int vii = qq.front();
        qq.pop();
//        cout << "Attr iteration, size: " << bigA->size() << " qq size: " << qq.size() << " with in edges: " << game->in_edges[vii].size() << " \n";

        for(auto & i : game->in_edges[vii]){
            int vi = target(i);
            if(bigV->find(vi) == bigV->end()) // vertex not in the playing area anymore
                continue;
            bool attracted = true;
            if(game->owner[vi] != player){
                for(auto & j : game->out_edges[vi]){
                    attracted = (attracted && (bigV->find(target(j)) == bigV->end()));
                }
            }
            if(!attracted)
                continue;
//            cout << vi << ",";
            bigA->insert(vi);
            bigV->erase(vi);
            qq.push(vi);
        }
    }
    cout << "\nAttr end, size: " << bigA->size() << "\n";
}
#else
void zlnk::attrQueue(int player, unordered_set<int> *bigA, vector<Subset> *ac) {
    cout << "Attr start, size: " << bigA->size() << "\n";
    queue<int> qq;
    for (const auto& vi : *bigA) {
        qq.push(vi);
    }
    while(!qq.empty())
    {
        int vii = qq.front();
        qq.pop();
//        cout << "Attr iteration, size: " << bigA->size() << " qq size: " << qq.size() << " with in edges: " << game->in_edges[vii].size() << " \n";

        for(int i = 0;i<game->in_edges[vii].size();i++){
            int vi = target(game->in_edges[vii][i]);
            int gi = guard_index(game->in_edges[vii][i]);
            if(bigV->find(vi) == bigV->end()) // vertex not in the playing area anymore
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
            if(bigA->find(vi) == bigA->end()){
                // add vertex to attracted set
                bigA->insert(vi);
                (*ac)[vi] = attracted;
            } else {
                (*ac)[vi] |= attracted;
            }
            (*vc)[vi] -= attracted;
            if((*vc)[vi] == emptyset)
            {
                bigV->erase(vi);
            }
            qq.push(vi);
        }
    }
    cout << "Attr end, size: " << bigA->size() << "\n";
}
#endif


void zlnk::test() {
    auto * bigA = new unordered_set<int>;
    vector<Subset> * ac = new vector<Subset>(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++)
        (*ac)[i] = emptyset;
    bigA->insert(12);
    (*ac)[12] = game->bigC;
    attr(1, bigA, ac);
}

#ifdef SINGLEMODE
/**
 * Assume empty winning sets (bigV empty and vc all 0's)
 * @param W0bigV
 * @param W0vc
 * @param W1bigV
 * @param W1vc
 */
void zlnk::solve(unordered_set<int> *W0bigV, vector<Subset> *W0vc, unordered_set<int> *W1bigV, vector<Subset> *W1vc) {
    if(bigV->empty()) return;

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

    unordered_set<int> *WMebigV;
    unordered_set<int> *WOpbigV;
    if(player == 0){
        WMebigV = W0bigV;
        WOpbigV = W1bigV;
    } else {
        WMebigV = W1bigV;
        WOpbigV = W0bigV;
    }



    auto * bigA = new unordered_set<int>;

    auto * subBigV = new unordered_set<int>;
    getVCWithPrio(bigA, nullptr, h);

    *subBigV = *bigV;
    zlnk subgame(game, subBigV, nullptr);
    subgame.attr(player, bigA, nullptr);
    cout << "\nDown1\n";
    fflush(stdout);
    subgame.solve(W0bigV, W0vc, W1bigV, W1vc);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    fflush(stdout);
    if(WOpbigV->empty()){
        unify(WMebigV, nullptr, bigA, nullptr);
        return;
    } else {
        // clone content and wipe winning sets
        *bigA = *WOpbigV;
        W0bigV->clear();
        W1bigV->clear();

        *subBigV = *bigV;
        zlnk subgame2(game, subBigV, vc);
        subgame2.attr(1-player, bigA, nullptr);
        cout << "\nDown2\n";
        fflush(stdout);
        subgame2.solve(W0bigV, W0vc, W1bigV, W1vc);
        attracting += subgame.attracting;
        cout << "\nUp2\n";
        fflush(stdout);
        unify(WOpbigV, nullptr, bigA, nullptr);
    }
}
#else
/**
 * Assume empty winning sets (bigV empty and vc all 0's)
 * @param W0bigV
 * @param W0vc
 * @param W1bigV
 * @param W1vc
 */
void zlnk::solve(unordered_set<int> *W0bigV, vector<Subset> *W0vc, unordered_set<int> *W1bigV, vector<Subset> *W1vc) {
    if(bigV->empty()) return;

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

    unordered_set<int> *WMebigV;
    vector<Subset> *WMevc;
    unordered_set<int> *WOpbigV;
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



    auto * bigA = new unordered_set<int>;
    vector<Subset> * ac = new vector<Subset>(game->n_nodes);

    auto * subBigV = new unordered_set<int>;
    vector<Subset> * subvc = new vector<Subset>(game->n_nodes);

    for(int i = 0;i<game->n_nodes;i++){
        (*ac)[i] = emptyset;
    }
    getVCWithPrio(bigA, ac, h);

    *subBigV = *bigV;
    *subvc = *vc;
    zlnk subgame(game, subBigV, subvc);
    subgame.attr(player, bigA, ac);
    cout << "\nDown1\n";
    fflush(stdout);
    subgame.solve(W0bigV, W0vc, W1bigV, W1vc);
    attracting += subgame.attracting;
    cout << "\nUp1\n";
    fflush(stdout);
    if(WOpbigV->empty()){
        unify(WMebigV, WMevc, bigA, ac);
        return;
    } else {
        // clone content and wipe winning sets
        *bigA = *WOpbigV;
        *ac = *WOpvc;
        W0bigV->clear();
        W1bigV->clear();
        for(int i = 0;i<game->n_nodes;i++){
            (*W0vc)[i] = emptyset;
            (*W1vc)[i] = emptyset;
        }

        *subBigV = *bigV;
        *subvc = *vc;
        zlnk subgame2(game, subBigV, vc);
        subgame2.attr(1-player, bigA, ac);
        cout << "\nDown2\n";
        fflush(stdout);
        subgame2.solve(W0bigV, W0vc, W1bigV, W1vc);
        attracting += subgame.attracting;
        cout << "\nUp2\n";
        fflush(stdout);
        unify(WOpbigV, WOpvc, bigA, ac);
    }
}
#endif

tuple<int, int> zlnk::getHighLowPrio() {
    int highest = 0;
    int lowest = INT32_MAX;
    for (const auto& vi : *bigV) {
        int prio = game->priority[vi];
        if(prio < lowest) lowest = prio;
        if(prio > highest) highest = prio;
    }
    return make_tuple(highest, lowest);
}

void zlnk::getVCWithPrio(unordered_set<int> *bigA, vector<Subset> *ac, int prio) {
    for (const auto& vi : *bigV) {
        if(game->priority[vi] == prio){
            bigA->insert(vi);
#ifndef SINGLEMODE
            (*ac)[vi] = (*vc)[vi];
#endif
        }
    }
}

void zlnk::unify(unordered_set<int> *bigA, vector<Subset> *ac, unordered_set<int> *bigB, vector<Subset> *bc) {
    for (const auto& vi : *bigB) {
        bigA->insert(vi);
#ifndef SINGLEMODE
        (*ac)[vi] |= (*bc)[vi];
#endif
    }
}

#ifdef SINGLEMODE
void zlnk::removeFromBigV(unordered_set<int> * bigA, vector<Subset> *ac){
    for (const auto& vi : *bigA) {
        bigV->erase(vi);
    }
}
void zlnk::removeFromBigV(int i, Subset c) {
    bigV->erase(i);
}

#else
void zlnk::removeFromBigV(unordered_set<int> * bigA, vector<Subset> *ac){
    for (const auto& vi : *bigA) {
        removeFromBigV(vi, (*ac)[vi]);
    }
}

void zlnk::removeFromBigV(int i, Subset c) {
    (*vc)[i] -= c;
    if((*vc)[i] == emptyset)
    {
        bigV->erase(i);
    }
}

#endif