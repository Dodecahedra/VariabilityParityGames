//
// Created by sjef on 5-6-19.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include "BDD/bddObj.h"
#include "Game.h"


void Game::set_n_nodes(int nodes) {
    n_nodes = nodes;
    out_edges = new std::vector<std::tuple<int,int>>[n_nodes];
    in_edges = new std::vector<std::tuple<int,int>>[n_nodes];
    priority.resize(n_nodes);
    owner.resize(n_nodes);
    declared.resize(n_nodes);
    for(int i = 0;i<n_nodes;i++)
        declared[i] = false;
}


Game::Game() {

}


void Game::parseGameFromFile(const string& filename) {
    int c = 0;
    std::ifstream infile(filename);

    while (infile.good())
    {
        char s[PARSER_LINE_SIZE];
        infile.getline(s, PARSER_LINE_SIZE,';');
        if(c == 0){
            // create bigC
            cout << "Found confs: " << s << '\n';
            parseConfs(s);
            c++;
        } else if(c == 1)
        {
            cout << "Found game: " << s << '\n';
            // init with parity
            parseInitialiser(s);
            c++;
        } else {
            // add vertex
            if(strlen(s) > 0)
                parseVertex(s);
        }
    }
    if (!infile.eof()) {
        throw std::string("could not open file");
    }
}

void Game::parseConfs(char * line) {
    while(*line == '\n' || *line == '\t' ||*line == ' ')
        line++;
    if(strncmp(line, "confs ",6) != 0) throw std::string("Expected confs");
    char conf[PARSER_LINE_SIZE];
    int i = 6;
    char c;
    do{
        c = line[i++];
    } while(c != '\0' && c != '+');
    bm_n_vars = i - 7;
    bm_vars.resize(bm_n_vars);
    bm = new bddMgr(bm_n_vars);
    for(int i = 0;i<bm_n_vars;i++) {
        bm_vars[i] = bm->bddVar(i);
    }
    parseConfSet(line, 6, &bigC);
}

void Game::parseInitialiser(char *line) {
    while(*line == '\n' || *line == '\t' ||*line == ' ')
        line++;
    if(strncmp(line, "parity ",7) != 0) throw std::string("Expected parity");
    int parity = atoi(line + 7);
    set_n_nodes(parity);
}


int Game::parseConfSet(const char *line, int i, BDD *result) {
    *result = bm->getZero();
    BDD entry = bm->getOne();
    int var = 0;
    char c;
    do
    {
        c = line[i++];
        if(c == '0'){
            if(var > bm_n_vars) throw std::string("Too many bits");
            entry = entry & ~bm_vars[var];
            var++;
        } else if(c == '1'){
            if(var > bm_n_vars) throw std::string("Too many bits");
            entry = entry & bm_vars[var];
            var++;
        } else if(c == '-'){
            if(var > bm_n_vars) throw std::string("Too many bits");
            var++;
        } else {
            *result = *result | entry;
            entry = bm->getOne();
            var = 0;
        }
    } while(c =='0' || c == '1' || c == '-' || c == '+');
    return i;
}

void Game::dumpSet(BDD *bdd, BDD t, char * p, int var) {
    if(var == bm_n_vars)
    {
        p[var]  = '\0';
        if(!((*bdd & t) == bm->getZero())){
            cout << p << '\n';
        }
    } else {
        BDD t1 = t & bm_vars[var];
        p[var] = '1';
        dumpSet(bdd, t1, p , var + 1);
        p[var] = '0';
        BDD t2 = t & ~bm_vars[var];
        dumpSet(bdd, t2, p, var + 1);
    }
}

void Game::parseVertex(char *line) {
    while(*line == '\n' || *line == '\t' ||*line == ' ')
        line++;
    int index;
    int i;
    i = readUntil(line, ' ');
    index = atoi(line);
    if(declared[index]) throw std::string("Already declared vertex " + std::to_string(index));
    line += i + 1;
    i = readUntil(line, ' ');
    priority[index] = atoi(line);
    line += i + 1;
    i = readUntil(line, ' ');
    owner[index] = atoi(line);
    line += i + 1;


//    cout << "Vertex with index: " << index << " and prio: " << priority[index] << " and owner: " << owner[index] << "\n";
    while(*line != '\0')
    {
        i = readUntil(line, '|');
        int target = atoi(line);
        line += i + 1;

        int guardindex = edge_guards.size();
        edge_guards.resize(guardindex + 1);
        i = parseConfSet(line, 0,&edge_guards[guardindex]);
        edge_guards[guardindex] = edge_guards[guardindex] & bigC;
        int outindex = out_edges[index].size();
        out_edges[index].resize(outindex + 1);
        out_edges[index][outindex] = std::make_tuple(target, guardindex);

        int inindex = in_edges[target].size();
        in_edges[target].resize(inindex+1);
        in_edges[target][inindex] = std::make_tuple(index, guardindex);
        //cout<< "with edge to " << target << " allowing: ";
        //dumpSet(&edge_guards[guardindex], bm->getOne(), new char[bm_n_vars+1], 0);
        line += i;
    }
    declared[index] = true;
}

int Game::readUntil(const char * line, char delim){
    int i = 0;
    while(*(line + i) != delim)
        i++;
    return i;
}
