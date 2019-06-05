//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_GAME_H
#define VPGSOLVER_GAME_H
#define PARSER_LINE_SIZE  16384
#include <vector>
#include <tuple>
#include "BDD/bddObj.h"

class Game {
public :
    bddMgr * bm;
    std::vector<BDD> bm_vars;
    int bm_n_vars;
    BDD bigC;
    int n_nodes;
    std::vector<int> priority;
    std::vector<int> owner;
    std::vector<bool> declared;
    std::vector<std::tuple<int,int>> *out_edges;
    std::vector<std::tuple<int,int>> *in_edges;
    std::vector<BDD> edge_guards;

    Game();
    void set_n_nodes(int nodes);

    void parseGameFromFile(const string& filename);
    void parseConfs(char * line);
    void parseInitialiser(char* line);
    int parseConfSet(const char * line, int i, BDD * result);
    void parseVertex(char * line);

    void dumpSet(BDD * bdd, BDD t, char * p, int var);
    int readUntil(const char * line, char delim);
};


#endif //VPGSOLVER_GAME_H
