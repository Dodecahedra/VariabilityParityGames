//
// Created by sjef on 11-7-19.
//

#include "MBR.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <chrono>
#include "FPIte.h"

vector<Subset> MBR::winningConf;
vector<VertexSet> MBR::winningVertices;

MBR::MBR(Game *game, Subset *conf, vector<bool> *edgeenabled, VertexSet *P0, VertexSet *VP1, int feature) {
    this->game = game;
    this->conf = conf;
    this->edgeenabled = edgeenabled;
    this->P0 = P0;
    this->VP1 = VP1;
    this->feature = feature;
}

MBR::MBR(Game *game) {
    this->game = game;
    this->conf = new Subset;
    *this->conf = game->bigC;

    this->P0 =  new VertexSet();
    this->P0->resize(game->n_nodes);
    this->VP1 =  new VertexSet();
    this->VP1->resize(game->n_nodes);
    fill(this->VP1->begin(), this->VP1->end(), true);
    this->edgeenabled = new vector<bool>(game->edge_guards.size());
    fill(this->edgeenabled->begin(), this->edgeenabled->end(), true);
    this->feature =  0;
}

void MBR::solve() {
    if(this->metric_output){
        if(this->measured == nullptr)
            this->measured = new bintree<vector<int>>();
        this->measured->value = new vector<int>(4);

#ifdef subsetbdd
        (*this->measured->value)[0] = std::count_if(P0->begin(), P0->end(), [](bool b){return b;});
        (*this->measured->value)[1] = game->n_nodes - std::count_if(VP1->begin(), VP1->end(), [](bool b){return b;});
#endif
    }
    if(this->feature == game->bm_n_vars){
        int i = winningConf.size();
        winningConf.resize(i+1);
        winningVertices.resize(i+1);
        winningConf[i] = *conf;
        winningVertices[i].resize(game->n_nodes);
        FPIte fpite(game, P0, VP1, edgeenabled, &(winningVertices[i]));
        fpite.solvelocal = solvelocal;

        if(this->metric_output)
        {
            auto start = std::chrono::high_resolution_clock::now();
            fpite.solve();
            auto end = std::chrono::high_resolution_clock::now();
            (*this->measured->value)[2] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            cout << "Assisted " << (*this->measured->value)[0] + (*this->measured->value)[1] << " with time " << (*this->measured->value)[2] << endl;
        } else {
            fpite.solve();
        }
        return;
    }

    vector<bool> pessimisticedges0(edgeenabled->size());
    vector<bool> pessimisticedges1(edgeenabled->size());

    createPessimisticGames(&pessimisticedges0, &pessimisticedges1);

    VertexSet W0(game->n_nodes);
    auto * P0b = new VertexSet;
    auto * VP1b = new VertexSet;
//
//    if(this->feature % 2 == 0) {
//    attr(0, P0, &pessimisticedges0);
//    attr(1, VP1, this->edgeenabled);
        auto *fpite0 = new FPIte(game, P0, VP1, &pessimisticedges0, &W0);
        if (this->metric_output) {
            auto start = std::chrono::high_resolution_clock::now();
            fpite0->solve();
            auto end = std::chrono::high_resolution_clock::now();
            (*this->measured->value)[2] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            cout << "Assisted " << (*this->measured->value)[0] + (*this->measured->value)[1] << " with time " << (*this->measured->value)[2] << endl;
        } else {
            fpite0->solve();
        }
        *P0 = W0;

//    *P0b = *P0;
//    attr(0, P0b, &pessimisticedges1);
//    auto * fpite1 = new FPIte(game, P0b, VP1, &pessimisticedges1, &W0);
        auto *fpite1 = new FPIte(game, P0, VP1, &pessimisticedges1, &W0);
        if (this->metric_output) {
            auto start = std::chrono::high_resolution_clock::now();
            fpite1->solve();
            auto end = std::chrono::high_resolution_clock::now();
            (*this->measured->value)[3] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            cout << "Assisted " << (*this->measured->value)[0] + (*this->measured->value)[1] << " with time " << (*this->measured->value)[3] << endl;
        } else {
            fpite1->solve();
        }
        *VP1 = W0;

        delete fpite0;
        delete fpite1;

        bool done;
        if(this->solvelocal){
            //Local:
            done = (*P0)[game->reindexedNew[0]] || !(*VP1)[game->reindexedNew[0]];
        } else {
            //Global:
            done = *P0 == *VP1;
        }
        if (done) {
            int i = winningConf.size();
            winningConf.resize(i + 1);
            winningVertices.resize(i + 1);
            winningConf[i] = *conf;
            winningVertices[i] = *P0;
            return;
        }
//    }

    auto * confa = new Subset;
    auto * confb = new Subset;
    do {
        *confa = *conf;
        parition(confa, confb);
        this->feature++;
        if(this->feature >  game->bm_n_vars)
            throw "Empty conf sets found in partitioning";
    } while (*confa == emptyset || *confb == emptyset);

    fill(pessimisticedges0.begin(), pessimisticedges0.end(), false);
    createSubGames(confa, &pessimisticedges0);
    createSubGames(confb, edgeenabled);

    *P0b = *P0;
    *VP1b = *VP1;

//    VertexSet dummy(P0->size());
//
//    fill(dummy.begin(), dummy.end(), false);
//    fill(dummy.begin(), dummy.end(), true);
//    fill(dummy.begin(), dummy.end(), false);
//    fill(dummy.begin(), dummy.end(), true);
//    fill(P0->begin(), P0->end(), false);
//    fill(P0b->begin(), P0b->end(), false);
//    fill(VP1->begin(), VP1->end(), true);
//    fill(VP1b->begin(), VP1b->end(), true);

    MBR ma(game, confa, &pessimisticedges0, P0, VP1, feature);
    MBR mb(game, confb, edgeenabled, P0b, VP1b, feature);
    if(this->metric_output){
        ma.metric_output = metric_output;
        mb.metric_output = metric_output;
        this->measured->left = new bintree<vector<int>>();
        this->measured->right = new bintree<vector<int>>();
        ma.measured = this->measured->left;
        mb.measured = this->measured->right;
    }
    ma.solvelocal = solvelocal;
    mb.solvelocal = solvelocal;
    ma.solve();
    mb.solve();

    delete P0b;
    delete VP1b;
}

void MBR::createPessimisticGames(vector<bool> *pessimisticedges0, vector<bool> *pessimisticedges1) {
    vector<bool> * pess;
    vector<bool> * opt;
    for(int i = 0;i<edgeenabled->size();i++){
        if(!(*edgeenabled)[i])
            continue;
        if(game->owner[game->edge_origins[i]] == 0){
            pess = pessimisticedges0;
            opt = pessimisticedges1;
        } else{
            pess = pessimisticedges1;
            opt = pessimisticedges0;
        }
        (*opt)[i] = true;
        Subset g = *conf;
        g -= game->edge_guards[i];
        (*pess)[i] = g == emptyset;
    }
}

void MBR::parition(Subset *org, Subset *part) {
    *part = *org;
    *org &= game->bm_vars[feature];
    *part -= game->bm_vars[feature];
}

void MBR::createSubGames(Subset *confP, vector<bool> *edgeenabledP) {
    for(int i = 0;i<edgeenabled->size();i++) {
        if (!(*edgeenabled)[i])
            continue;
        Subset g = *confP;
        g &= game->edge_guards[i];
        (*edgeenabledP)[i] = !(g == emptyset);
    }
}

void MBR::attr(int player, VertexSet *U, vector<bool> * edgeenabledvector) {
    vector<bool> countinitialized;
    vector<int> countincoming;
    countincoming.resize(game->n_nodes);
    countinitialized.resize(game->n_nodes);
    queue<int> qq;
    for(int i = 0;i<game->n_nodes;i++)
        if((*U)[i])
            qq.push(i);
    while(!qq.empty())
    {
        int vii = qq.front();
        qq.pop();

        for(auto & i : game->in_edges[game->reindexedOrg[vii]]){
            if(!(*edgeenabledvector)[guard_index(i)])
                continue;
            int vi = game->reindexedNew[target(i)];
            if((*U)[vi]) // vertex already attracted
                continue;
            bool attracted = true;
            if(game->owner[game->reindexedOrg[vi]] != player){
                if(!countinitialized[vi]){
                    countinitialized[vi] = true;
                    countincoming[vi] = 0;
                    for(auto & j : game->out_edges[game->reindexedOrg[vi]])
                        if((*edgeenabledvector)[guard_index(j)])
                            countincoming[vi]++;
                }
                if((--countincoming[vi]) > 0){
                    attracted = false;
                }
//                attracted = false;
            }
            if(!attracted)
                continue;
            (*U)[vi] = true;
            qq.push(vi);
        }
    }
}

void MBR::printMeasurements(ostream * output) {
    *output << "digraph timetree {" << endl;
    int nodes = printNode(output, this->measured,0);
    *output << '}' << endl;
    cout << "nr of nodes: " << nodes << endl;
    cout << "Solved " << ((nodes+1)*3/2-2) << " games" << endl;
}

int MBR::printNode(ostream *output, bintree<vector<int>> *node, int c) {
    if(node == nullptr)
        return -1;
    *output << 'n' << c++ << " [label=\"";
    *output << "V0: " << (*node->value)[0] << '/' << this->game->n_nodes << endl;
    *output << "V1: " << (*node->value)[1] << '/' << this->game->n_nodes << endl;
    *output << "Solve 1: " << (*node->value)[2] << "ns" << endl;
    *output << "Solve 2: " << (*node->value)[3] << "ns";
    *output << "\"];" << endl;
    int l = printNode(output, node->left, c);
    if(l == -1) {
        int r = printNode(output, node->right, c);
        if(r != -1){
            *output << 'n' << (c-1) << " -> n" << c << ';' << endl;
            return r;
        }
        return c;
    } else {
        int r = printNode(output, node->right, l);
        *output << 'n' << (c-1) << " -> n" << c << ';' << endl;
        if(r != -1){
            *output << 'n' << (c-1) << " -> n" << l << ';' << endl;
            return r;
        }
        return l;
    }
}
