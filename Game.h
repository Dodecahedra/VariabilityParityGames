//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_GAME_H
#define VPGSOLVER_GAME_H
#define PARSER_LINE_SIZE  16384
#include <vector>
#include <tuple>
#include <unordered_set>
#include "bdd.h"

using namespace std;

#define target(a) std::get<0>(a)
#define guard_index(a) std::get<1>(a)

class Game {
public :
    std::vector<bdd> bm_vars;
    int bm_n_vars;
    bdd bigC;
    int n_nodes;
    std::vector<int> priority;
    std::vector<int> owner;
    std::vector<bool> declared;
    std::vector<std::tuple<int,int>> *out_edges;
    std::vector<std::tuple<int,int>> *in_edges;
    std::vector<bdd> edge_guards;

    Game();
    void set_n_nodes(int nodes);

    void parseGameFromFile(const string& filename);
    void parseConfs(char * line);
    void parseInitialiser(char* line);
    int parseConfSet(const char * line, int i, bdd * result);
    void parseVertex(char * line);

    void dumpSet(bdd * dumpee, bdd t, char * p, int var);
    void printCV(unordered_set<int> *bigV, vector<bdd> *vc);
    void printCV(unordered_set<int> *bigV, vector<bdd> *vc, bdd t, char * p, int var);
    int readUntil(const char * line, char delim);
};


#endif //VPGSOLVER_GAME_H
