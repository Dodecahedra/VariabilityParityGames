//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_GAME_H
#define VPGSOLVER_GAME_H
#define PARSER_LINE_SIZE  16777216
#include <vector>
#include <tuple>
#include <unordered_set>
#include <map>
#include <iostream>
#include "conf.h"

using namespace std;

#define target(a) std::get<0>(a)
#define guard_index(a) std::get<1>(a)

class Game {
public :
    std::vector<Subset> bm_vars;
    int bm_n_vars;
    Subset bigC;
    int n_nodes;
    vector<int> priority;
    vector<unordered_set<int>> priorityI;
    vector<int> owner;
    vector<bool> declared;
    vector<std::tuple<int,int>> *out_edges;
    vector<std::tuple<int,int>> *in_edges;
    vector<Subset> edge_guards;
    vector<int> edge_origins;
    vector<int> reindexedNew;
    vector<int> reindexedOrg;
    vector<int> reindexPCutoff;
    vector<vector<int>> orgvertices;
    int winningfor0 = -1;

    bool parsePG = false;
    bool compressvertices = false;

    bool specificvarlast = false;
    int specificvar;

    map<string, Subset> parseCache;

    Game();
    void set_n_nodes(int nodes);

    void parseVPGFromFile(const string &filename);
    void parseVPGFromFile(const string &filename, const char *specificconf);
    void parseConfs(char * line);
    void parseInitialiser(char* line);
    int parseConfSet(const char * line, int i, Subset * result);
    void parseVertex(char * line);

    void dumpSet(Subset * dumpee, Subset t, char * p, int var);
    void printCV(unordered_set<int> *bigV, vector<Subset> *vc, bool fulloutut);
    void printCV(unordered_set<int> *bigV, vector<Subset> *vc, Subset t, char * p, int var, bool fulloutput);
    int readUntil(const char * line, char delim);

    void buildInEdges();
    void compressPriorities();
    void compressVertices();
    void moveVertexInto(int v1, int v2, bool v1Ev2);
    void movePriorities(int from, int to);
    void reindexVertices();
    int findVertexWinningForVertex(int v);
    int findVertexWinningFor0();

    void parsePGFromFile(const string &filename);

    void writePG(ostream * output);
    void writePG(ostream * output, Subset conf);
    void findAllElements(Subset s, vector<tuple<Subset, string>> * result);
    void findAllElements(Subset s, vector<tuple<Subset, string>> * result, char * p, int var);
};


#endif //VPGSOLVER_GAME_H
