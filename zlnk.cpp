//
// Created by sjef on 5-6-19.
//

#include "zlnk.h"

zlnk::zlnk(Game *game) {
    unordered_set<int> * bigV = new unordered_set<int>(game->n_nodes);
    vector<BDD> * vc = new vector<BDD>(game->n_nodes);
    for(int i = 0;i<game->n_nodes;i++){
        (*vc)[i] = game->bigC;
    }
    zlnk(game, bigV, vc);
}
zlnk::zlnk(Game * game, unordered_set<int> * bigV, vector<BDD> * vc){
    this->game = game;
    this->bigV = bigV;
    this->vc = vc;
}

void zlnk::attr(int player, std::vector<tuple<int, BDD>> a, int ia) {
    bool workdone = true;
    while(workdone){
        for(int i = 0;i<ia;i++){
            int v = std::get<0>(a[i]);
            BDD vc = std::get<1>(a[i]);
            for(int ii = 0;ii < game->in_edges[i].size();i++){
                int f = std::get<0>(game->in_edges[i][ii]);
                int eci = std::get<1>(game->in_edges[i][ii]);
                BDD ec = game->edge_guards[eci];
                if(game->owner[f] == player){
                    BDD vcn = ec & vc; // & fc
                }
            }
        }
    }
}
