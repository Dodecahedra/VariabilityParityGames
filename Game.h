/***********************************************************************************************************************
 * In this project only algorithms and datastructures are implemented properly.
 *
 * The code in this file is not optimized and not up to standards; it is sufficient only for experimental applications
 * but not for any real application.
 **********************************************************************************************************************/
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
#define edge_index(a) std::get<1>(a)

/**
 * Represent VPG using a double edge relation. Every edge has an edge index.
 * Every edge index is mapped to a set of configurations guarding the edge.
 */
class Game {
public :
    std::vector<ConfSet> bm_vars;
    int bm_n_vars;
    ConfSet bigC;
    int n_nodes;
    vector<int> priority;
    vector<unordered_set<int>> priorityI;
    vector<int> owner;
    vector<bool> declared;
    vector<std::tuple<int,int>> *out_edges;
    vector<std::tuple<int,int>> *in_edges;
    vector<ConfSet> edge_guards;
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

    map<string, ConfSet> parseCache;

    Game();
    void set_n_nodes(int nodes);

    void parseVPGFromFile(const string &filename);
    void parseVPGFromFile(const string &filename, const char *specificconf);
    void parseConfs(char * line);
    void parseInitialiser(char* line);
    int parseConfSet(const char * line, int i, ConfSet * result);
    void parseVertex(char * line);

    void dumpSet(ConfSet * dumpee, ConfSet t, char * p, int var);
    void printCV(VertexSetZlnk *bigV, vector<ConfSet> *vc, bool fulloutut);
    void printCV(VertexSetZlnk *bigV, vector<ConfSet> *vc, ConfSet t, char * p, int var, bool fulloutput);
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
    void writePG(ostream * output, ConfSet conf);
    void findAllElements(ConfSet s, vector<tuple<ConfSet, string>> * result);
    void findAllElements(ConfSet s, vector<tuple<ConfSet, string>> * result, char * p, int var);
};


#endif //VPGSOLVER_GAME_H
