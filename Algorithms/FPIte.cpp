
//
// Created by sjef on 11-7-19.
//

#include <random>
#include "FPIte.h"
#include <algorithm>
#include <iostream>
#include <chrono>

#define targetIsIn(t) ZZ[t]

FPIte::FPIte(Game *game, VertexSetFPIte *P0, VertexSetFPIte *VP1, VertexSetFPIte * W0) {
    this->P0 = P0;
    this->VP1 = VP1;
    this->W0 = W0;

    this->game = game;
}

FPIte::FPIte(Game *game) {
    this->game = game;
    this->P0 =  new VertexSetFPIte();
    this->P0->resize(game->n_nodes);
    this->VP1 =  new VertexSetFPIte();
    this->VP1->resize(game->n_nodes);
    fill(this->VP1->begin(), this->VP1->end(), true);
    this->W0 = new VertexSetFPIte();
    this->W0->resize(game->n_nodes);
}

void FPIte::init(int i, int ie) {
    if(i >= d || ie >= d || i < 0 || ie < 0)
        return;
    if(i % 2 == 0)
        copyWithPrio(&ZZ, VP1, i, ie);
    else
        copyWithPrio(&ZZ, P0, i, ie);
}

void FPIte::diamondbox(VertexSetFPIte *Z, int maxprio) {
    dbs_executed++;
    for(int p = 0;p <= maxprio;p++) {
        for (int v = game->reindexPCutoff[p]; v < game->reindexPCutoff[p + 2]; v++) {
            if (targetIsIn(v) != targetWasIn[v]) {
                for (auto &edge : game->in_edges[v]) {
                    //predesseccor of v
                    int t = target(edge);
                    // Vertices in P0 are always in and vertices in P1 are always out, so we don't need to reconsider them
                    if ((*P0)[t] || !(*VP1)[t])
                        continue;
                    verticesconsidered++;

                    if (game->owner[t] == 0) { //diamond
                        if (targetIsIn(v))
                            // v was not ZZ but now is, so we increase the number of successors that are in ZZ by 1
                            edgecount[t]++;
                        else
                            // v was in ZZ but no longer is, so we decrease the number of successors that are in ZZ by 1
                            edgecount[t]--;
                        // as long as there is 1 successor in ZZ does the diamond formula hold for t
                        (*Z)[t] = edgecount[t] > 0;
                    } else { //box
                        if (targetIsIn(v))
                            // v was not in ZZ but now is, so we decrease the number of successors not in ZZ by 1
                            edgecount[t]--;
                        else
                            // v was in ZZ but no loner is, so we increase the number of successor ot in ZZ by 1
                            edgecount[t]++;
                        // only if all successors are in ZZ does the box formula hold for t
                        (*Z)[t] = edgecount[t] == 0;
                    }
                }
                targetWasIn[v] = targetIsIn(v);
            }
        }
    }
}

void FPIte::diamondbox(VertexSetFPIte *Z) {
    edgecount.resize(game->n_nodes);
    targetWasIn.resize(game->n_nodes);
    for(int v = 0;v<game->n_nodes;v++){
        bool in;
        // Vertices in P0 and outside VP1 are not considered
        if((*P0)[v]) {
            in = true;
            targetWasIn[v] = true;
        } else if(!(*VP1)[v])
            in = false;
        else {
            targetWasIn[v] = targetIsIn(v);
            if (game->owner[v] == 0) {//diamond
                in = false;
                // do an initial count of the number of successors that are in ZZ
                for(auto edge : game->out_edges[v]){
                    int t = target(edge);
                    if (targetIsIn(t)) {
                        in = true; // if at least one successor is in then the diamond formula holds for v
                        edgecount[v]++;
                    }
                }
            } else {//box
                in = true;
                // do an initial count of the number of successors that are not in ZZ
                for(auto edge : game->out_edges[v]){
                    int t = target(edge);
                    if (!targetIsIn(t)) {
                        in = false; // if one successor is not in ZZ then the box formula does not hold for v
                        edgecount[v]++;
                    }
                }
            }
        }
        (*Z)[v] = in;
    }
}

void FPIte::solve() {
    d = game->priorityI.size();
    ZZ.resize(game->n_nodes);

    // initialize everything
    if((d-1) % 2 == 0){
        init(0,d - 1);
        init(1,d - 2);
    } else {
        init(1,d - 1);
        init(0,d - 2);
    }

    int i = 0;
    bool equal;
    bool first = true;
    do{
        if(solvelocal && i == d){
            // if we can solve local and we have just modified the left most fixed-point variable then maybe can terminate
            if((d - 1) % 2 == 0){
                // left most fixed-point variable is a greatest fixed-point
                if(!(*W0)[game->reindexedNew[game->findVertexWinningFor0()]])
                    break; // the vertex is no longer in W0 and therefore never will be so we can terminate
            } else {
                // left most fixed-point variable is a least fixed-point
                if((*W0)[game->reindexedNew[game->findVertexWinningFor0()]])
                    break; // the vertex is in W0 and therefore always will be so we can terminate
            }
        }
        if(first) {
            // initialize the edge count once
            diamondbox(W0);
            first = false;
        } else {
            diamondbox(W0, i - 1);
        }
        i = 1;
        if(!compareWithPrio(&ZZ, W0,0)){
            // keep calculating iteration variable 0
            copyWithPrio(&ZZ, W0, 0);
            equal = false;
        } else {
            // we have found a fixed-point for iteration variable 0.
            // we copy iteration variable 0 to 1 and if there is no difference between them copy 0 to 2 etc..
            equal = true;
            while (equal && i < d) {
                equal = compareWithPrio(&ZZ, W0, i);
                copyWithPrio(&ZZ, W0, i);
                i++;
            }
            // iteration variable 0 is copied into iteration variables 1 .. i-1
            // initialize the variables below i-1 and only those with parity different from i-1
            if(i % 2 == 0){
                init(0,i-2);
            } else {
                init(1,i-2);
            }
        }
        // Check if we were able to copy the iteration variable in variable d-1 and if those were also equal
    }while(!(i == d && equal));
}

FPIte::~FPIte() = default;

void FPIte::copyWithPrio(VertexSetFPIte *Z, VertexSetFPIte *ZP, int p) {
    int start = game->reindexPCutoff[p];
    Z->copy_n(ZP, start, game->reindexPCutoff[p+2] - start);
}

bool FPIte::compareWithPrio(VertexSetFPIte *Z, VertexSetFPIte *ZP, int p) {
    int start = game->reindexPCutoff[p];
    return Z->compare_n(ZP,start, game->reindexPCutoff[p+2] - start);
}

void FPIte::copyWithPrio(VertexSetFPIte *Z, VertexSetFPIte *ZP, int sp, int ep) {
    if(ep < sp)
        return;
    int start = game->reindexPCutoff[sp];
    Z->copy_n(ZP, start, game->reindexPCutoff[ep+2] - start);
}

void FPIte::setP0(char *P0string) {
    char * pEnd;
    do{
        int v = strtol(P0string, &pEnd, 10);
        (*P0)[v] = true;
        P0string = pEnd+1;
    } while (*P0string != '\0');
}

void FPIte::P0IsFull() {
    *VP1 = *P0;
}

void FPIte::unassist(int n) {
    vector<int> r (this->game->n_nodes);
    for(int i = 0;i<r.size();i++)
        r[i] = i;
    shuffle(r.begin(), r.end(), std::mt19937(std::random_device()()));
    for(int i = 0;i<n;i++)
        if((*P0)[i])
            (*P0)[i] = false;
        else
            (*VP1)[i] = true;
}
