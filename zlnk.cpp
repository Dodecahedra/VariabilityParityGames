//
// Created by sjef on 5-6-19.
//

#include <queue>
#include "zlnk.h"

zlnk::zlnk(Game *game) {
    unordered_set<int> * bigV = new unordered_set<int>(game->n_nodes);
    vector<bdd> * vc = new vector<bdd>(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++){
        bigV->insert(i);
        (*vc)[i] = game->bigC;
    }
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

zlnk::zlnk(Game * game, unordered_set<int> * bigV, vector<bdd> * vc){
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

void zlnk::attr(int player, unordered_set<int> *bigA, vector<bdd> *ac) {
    removeFromBigV(bigA, ac);
    attrQueue(player, bigA, ac);
}

/**
 *
 * @param player
 * @param bigA
 * @param ac It must hold that (*ac)[i] = bdd0 if i not in bigA
 */
void zlnk::attrNaive(int player, unordered_set<int> * bigA, vector<bdd> * ac){
    bool workdone = true;
    cout << "Attr start, size: " << bigA->size() << "\n";
    fflush(stdout);
    while(workdone){
        workdone = false;
        cout << "Attr iteration, size: " << bigA->size() << "\n";
        fflush(stdout);
        for (const auto& vi : *bigV) {
            bdd attracted;
            if(game->owner[vi] == player){
                attracted = bddfalse;
                for(int i = 0;i<game->out_edges[vi].size();i++){
                    int target = target(game->out_edges[vi][i]);
                    bdd guard = game->edge_guards[guard_index(game->out_edges[vi][i])];
                    if(guard == bddfalse || bigA->find(target) == bigA->end())
                        continue;
                    attracted = attracted | ((*vc)[vi] & (*ac)[target] & guard);
                }
            } else {
                attracted = (*vc)[vi];
                for(int i = 0;i<game->out_edges[vi].size();i++){
                    int target = target(game->out_edges[vi][i]);
                    bdd guard = game->edge_guards[guard_index(game->out_edges[vi][i])];
                    if(guard == bddfalse)
                        continue;
                    attracted = attracted & ((*ac)[target] | !guard);
                }
            }
            if(attracted == bddfalse)
                continue;
            workdone = true;
            if(bigA->find(vi) == bigA->end()){
                // add vertex to attracted set
                bigA->insert(vi);
                (*ac)[vi] = attracted;
            } else {
                (*ac)[vi] = (*ac)[vi] | attracted;
            }
            (*vc)[vi] = (*vc)[vi] & !attracted;
            if((*vc)[vi] == bddfalse)
            {
                bigV->erase(vi);
            }
        }
    }
}
void zlnk::attrQueue(int player, unordered_set<int> *bigA, vector<bdd> *ac) {
//    cout << "Attr start, size: " << bigA->size() << "\n";
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
            bdd guard = game->edge_guards[gi];
            bdd attracted;
            if(game->owner[vi] == player){
                attracted = (*vc)[vi] & (*ac)[vii] & guard;
            } else {
                attracted = (*vc)[vi];
                for(int i = 0;i<game->out_edges[vi].size();i++){
                    int target = target(game->out_edges[vi][i]);
                    bdd guard = game->edge_guards[guard_index(game->out_edges[vi][i])];
                    attracted = attracted & ((*ac)[target] | !guard);
                }
            }
            if(attracted == bddfalse)
                continue;
            if(bigA->find(vi) == bigA->end()){
                // add vertex to attracted set
                bigA->insert(vi);
                (*ac)[vi] = attracted;
            } else {
                (*ac)[vi] = (*ac)[vi] | attracted;
            }
            (*vc)[vi] = (*vc)[vi] & !attracted;
            if((*vc)[vi] == bddfalse)
            {
                bigV->erase(vi);
            }
            qq.push(vi);
        }
    }
}


void zlnk::test() {
    auto * bigA = new unordered_set<int>;
    vector<bdd> * ac = new vector<bdd>(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++)
        (*ac)[i] = bddfalse;
    bigA->insert(12);
    (*ac)[12] = game->bigC;
    attr(1, bigA, ac);
}

/**
 * Assume empty winning sets (bigV empty and vc all 0's)
 * @param W0bigV
 * @param W0vc
 * @param W1bigV
 * @param W1vc
 */
void zlnk::solve(unordered_set<int> *W0bigV, vector<bdd> *W0vc, unordered_set<int> *W1bigV, vector<bdd> *W1vc) {
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
    vector<bdd> *WMevc;
    unordered_set<int> *WOpbigV;
    vector<bdd> *WOpvc;
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
    vector<bdd> * ac = new vector<bdd>(game->n_nodes);

    auto * subBigV = new unordered_set<int>;
    vector<bdd> * subvc = new vector<bdd>(game->n_nodes);

    for(int i = 0;i<game->n_nodes;i++){
        (*ac)[i] = bddfalse;
    }
    getVCWithPrio(bigA, ac, h);

    *subBigV = *bigV;
    *subvc = *vc;
    zlnk subgame(game, subBigV, subvc);
    subgame.attr(player, bigA, ac);
    cout << "\nDown1\n";
    fflush(stdout);
    subgame.solve(W0bigV, W0vc, W1bigV, W1vc);
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
            (*W0vc)[i] = bddfalse;
            (*W1vc)[i] = bddfalse;
        }

        *subBigV = *bigV;
        *subvc = *vc;
        zlnk subgame2(game, subBigV, vc);
        subgame2.attr(1-player, bigA, ac);
        cout << "\nDown2\n";
        fflush(stdout);
        subgame2.solve(W0bigV, W0vc, W1bigV, W1vc);
        cout << "\nUp2\n";
        fflush(stdout);
        unify(WOpbigV, WOpvc, bigA, ac);
    }
}

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

void zlnk::getVCWithPrio(unordered_set<int> *bigA, vector<bdd> *ac, int prio) {
    for (const auto& vi : *bigV) {
        if(game->priority[vi] == prio){
            bigA->insert(vi);
            (*ac)[vi] = (*vc)[vi];
        }
    }
}

void zlnk::unify(unordered_set<int> *bigA, vector<bdd> *ac, unordered_set<int> *bigB, vector<bdd> *bc) {
    for (const auto& vi : *bigB) {
        bigA->insert(vi);
        (*ac)[vi] = (*ac)[vi] | (*bc)[vi];
    }
}

void zlnk::removeFromBigV(unordered_set<int> * bigA, vector<bdd> *ac){
    for (const auto& vi : *bigA) {
        removeFromBigV(vi, (*ac)[vi]);
    }
}

void zlnk::removeFromBigV(int i, bdd c) {
    (*vc)[i] = (*vc)[i] & !c;
    if((*vc)[i] == bddfalse)
    {
        bigV->erase(i);
    }
}
