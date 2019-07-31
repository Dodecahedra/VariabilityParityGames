//
// Created by sjef on 5-6-19.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <unordered_set>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::high_resolution_clock
#include <algorithm>    // std::shuffle

#include "Game.h"

void Game::set_n_nodes(int nodes) {
    n_nodes = nodes;
    out_edges = new std::vector<std::tuple<int,int>>[n_nodes];
    in_edges = new std::vector<std::tuple<int,int>>[n_nodes];
    priority.resize(n_nodes);
    owner.resize(n_nodes);
    declared.resize(n_nodes);


    reindexedNew.resize(this->n_nodes);
    reindexedOrg.resize(this->n_nodes);
    for(int v = 0;v < n_nodes;v++){
        reindexedOrg[v] = v;
        reindexedNew[v] = v;
    }
}


Game::Game() {

}


void Game::parseVPGFromFile(const string &filename, const char *specificconf) {
    int c = 0;
    if(parsePG){
        c++;
        parseConfs("confs 1;");
    }
    std::ifstream infile(filename);

    while (infile.good())
    {
        char s[PARSER_LINE_SIZE];
        infile.getline(s, PARSER_LINE_SIZE,';');
        if(c == 0){
            // create bigC
            cout << "Found confs: " << s << '\n';
            parseConfs(s);
            if(strlen(specificconf) > 0)
                parseConfSet(specificconf,0,&bigC);
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
    buildInEdges();
}

void Game::parseVPGFromFile(const string &filename) {
    parseVPGFromFile(filename, "");
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
    } while(c != '\0' && c != '+' && c != ';');
    bm_n_vars = i - 7;
    bm_vars.resize(bm_n_vars);
#ifdef subsetbdd
    bdd_init(2000000,2);
//    bdd_setcacheratio(200);
    bdd_setvarnum(bm_n_vars);
    vector<int> order;
    order.resize(bm_n_vars);
    for(i = 0;i<bm_n_vars;i++){
        order[i] = i;
    }
#ifdef randombddorder
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    shuffle (order.begin(), order.end(), default_random_engine(seed));
#else
    if(specificvarlast){
        order[bm_n_vars-1] = specificvar;
        order[specificvar] = bm_n_vars-1;
    }
#endif
    cout << "Bdd order: " ;
    for(i = 0;i<bm_n_vars;i++){
        cout << '[' << i << "]=" << order[i] << ", ";
    }
    cout << "\n";
    for(i = 0;i<bm_n_vars;i++) {
        bm_vars[order[i]] = bdd_ithvar(i);
    }
#endif
#ifdef subsetexplicit
    SubsetExplicit::size = (1 << (bm_n_vars));
    fullset.items.resize(SubsetExplicit::size);
    emptyset.items.resize(SubsetExplicit::size);
    for(int i =0;i<1<<bm_n_vars;i++) {
        fullset.items[i] = true;
        emptyset.items[i] = false;
    }
    for(int i = 0;i<bm_n_vars;i++) {
        bm_vars[i] = SubsetExplicit(bm_n_vars-i-1);
    }
#endif
    parseConfSet(line, 6, &bigC);
}

void Game::parseInitialiser(char *line) {
    while(*line == '\n' || *line == '\t' ||*line == ' ')
        line++;
    if(strncmp(line, "parity ",7) != 0) throw std::string("Expected parity");
    int parity = atoi(line + 7);
    set_n_nodes(parity);
}


int Game::parseConfSet(const char *line, int i, Subset *result) {
    if(parsePG){
        *result = fullset;
        return i+1;
    }
    *result = emptyset;
    Subset entry = fullset;
    int var = 0;
    char c;
    string * seq;
    do
    {
        if(var == 0){
            seq = new string(line+i);
            *seq = seq->substr(0, bm_n_vars);
            auto cacheite = parseCache.find(*seq);
            if(cacheite != parseCache.end())
            {
                *result |= cacheite->second;
                i += bm_n_vars;
                c = line[i++];
                continue;
            }
        }
        c = line[i++];
        if(c == 'F'){
            if(var != 0) throw std::string("Unexpected F");
            entry = emptyset;
        } else if (c == '0') {
            if (var > bm_n_vars) throw std::string("Too many bits");
            entry -= bm_vars[var];
            var++;
        } else if (c == '1') {
            if (var > bm_n_vars) throw std::string("Too many bits");
            entry &= bm_vars[var];
            var++;
        } else if (c == '-') {
            if (var > bm_n_vars) throw std::string("Too many bits");
            var++;
        } else {
            Subset * s = new Subset();
            *s = entry;
            parseCache.insert(pair<string, Subset>(*seq, *s));
            *result |= entry;
            entry = fullset;
            var = 0;
        }
    } while(c =='0' || c == '1' || c == '-' || c == '+' || c == 'F');
    return i;
}

void Game::dumpSet(Subset * dumpee, Subset t, char * p, int var) {
    if(var == bm_n_vars)
    {
        p[var]  = '\0';
        Subset result = *dumpee;
        result &= t;
        if(!(result == emptyset)){
            cout << p << ',';
        }
    } else {
        Subset t1 = t;
        t1 &= bm_vars[var];
        p[var] = '1';
        dumpSet(dumpee, t1, p , var + 1);
        p[var] = '0';
        Subset t2 = t;
        t2 -= bm_vars[var];
        dumpSet(dumpee, t2, p, var + 1);
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
    int p = atoi(line);
    priority[index] = p;
    if(p+1 > priorityI.size())
        priorityI.resize(p+1);
    priorityI[p].insert(index);
    line += i + 1;
    i = readUntil(line, ' ');
    owner[index] = atoi(line);
    line += i + 1;


    cout << "\nVertex with index: " << index << " and prio: " << priority[index] << " and owner: " << owner[index] << "\n";
    while(*line != '\0')
    {
        if(*line == ',')
            line++;
        int target;
        if(parsePG){
            i = readUntil(line, ',');
            target = atoi(line);
            line += i;
        } else {
            i = readUntil(line, '|');
            target = atoi(line);
            line += i + 1;
        }

        int guardindex = edge_guards.size();
        edge_guards.resize(guardindex + 1);
        edge_origins.resize(guardindex + 1);
        i = parseConfSet(line, 0,&edge_guards[guardindex]);
        edge_guards[guardindex] &= bigC;
        edge_origins[guardindex] = index;
        if(!(edge_guards[guardindex] == emptyset)){
            int outindex = out_edges[index].size();
            out_edges[index].resize(outindex + 1);
            out_edges[index][outindex] = std::make_tuple(target, guardindex);

//            int inindex = in_edges[target].size();
//            in_edges[target].resize(inindex+1);
//            in_edges[target][inindex] = std::make_tuple(index, guardindex);
//            cout<< "with edge to " << target << " allowing: ";
//            dumpSet(&edge_guards[guardindex], fullset, new char[bm_n_vars+1], 0);
        }
        line += i-1;
    }
    declared[index] = true;
}

int Game::readUntil(const char * line, char delim){
    int i = 0;
    while(*(line + i) != delim && *(line + i) != '\0')
        i++;
    return i;
}

void Game::printCV(unordered_set<int> *bigV, vector<Subset> *vc, Subset t, char * p, int var, bool fulloutput) {
    if(t == emptyset) return;
    if(var == bm_n_vars)
    {
        p[var]  = '\0';
        cout << "For product " << p << " the following vertices are in: ";
        if(fulloutput){
            for(const int& vi : *bigV) {
                Subset result = (*vc)[vi];
                result &= t;
                if(!((result) == emptyset)){
                    if(compressvertices){
                        for(auto j : orgvertices[vi])
                            cout << j << ',';
                    } else {
                        cout << vi << ',';
                    }
                }
            }
        } else {
            int vi = 0;
            if(compressvertices){
                bool found = false;
                for(int i = 0;i< n_nodes && !found;i++){
                    for(int j = 0;j < orgvertices[i].size() && !found;j++){
                        if(orgvertices[i][j] == 0){
                            vi = i;
                            found = true;
                        }
                    }
                }
            }
            Subset result = (*vc)[vi];
            result &= t;
            if(!((result) == emptyset)){
                cout << vi << ',';
            }
        }
        cout <<"\n";
        fflush(stdout);

    } else {
        Subset t1 = t;
        t1 &= bm_vars[var];
        p[var] = '1';
        printCV(bigV,vc, t1, p , var + 1, fulloutput);
        p[var] = '0';
        Subset t2 = t;
        t2 -= bm_vars[var];
        printCV(bigV, vc, t2, p, var + 1, fulloutput);
    }
}


void Game::printCV(unordered_set<int> *bigV, vector<Subset> *vc, bool fulloutput) {
#ifdef SINGLEMODE
    cout << "The following vertices are in: ";
    if(fulloutput){
        for(const int& vi : *bigV) {
            cout << vi << ',';
        }
    } else {
        if(bigV->find(0) != bigV->end())
            cout << "0,";
    }
    cout << "\n";
    fflush(stdout);
#else
    printCV(bigV, vc, bigC, new char[bm_n_vars+1], 0, fulloutput);
#endif
}

void Game::compressPriorities() {
    int pp = -1;
    for(int i = 0;i<priorityI.size();i++){
        if(!priorityI[i].empty())
        {
            if(pp == -1){
                pp = i;
                continue;
            }
            if(i - pp == 1){
                pp = i;
                continue;
            }
            if((i - pp) % 2 == 1){
                // odd
                movePriorities(i, ++pp);
                cout << "Move prio from " << i << " to " << pp << "\n";
            } else {
                // even
                movePriorities(i, pp);
                cout << "Move prio from " << i << " to " << pp << "\n";
            }
        }
    }
    priorityI.resize(pp+1);
}

void Game::movePriorities(int from, int to) {
    for(auto &e : priorityI[from]){
        priority[e] = to;
    }
    if(priorityI[to].empty()){
        priorityI[to] = priorityI[from];
    } else {
        priorityI[to].insert(priorityI[from].begin(), priorityI[from].end());
    }
    priorityI[from].clear();
}

void Game::reindexVertices() {
    reindexPCutoff.resize(this->priorityI.size() + 2);

    int c = 0;
    reindexPCutoff[0] = c;
    int i;
    for(i = 0;i < this->priorityI.size();i+=2){
        for(auto& v : priorityI[i]){
            reindexedNew[v] = c;
            reindexedOrg[c++] = v;
        }
        reindexPCutoff[i+2] = c;
    }
    reindexPCutoff[1] = c;
    for(i = 1;i < this->priorityI.size();i+=2){
        for(auto& v : priorityI[i]){
            reindexedNew[v] = c;
            reindexedOrg[c++] = v;
        }
        reindexPCutoff[i+2] = c;
    }
    vector<int> priority_old = priority;
    vector<int> owner_old  = owner;
    auto *out_edges_old = new vector<std::tuple<int,int>>[n_nodes];
    auto *in_edges_old = new vector<std::tuple<int,int>>[n_nodes];
    std::copy(out_edges, out_edges+n_nodes, out_edges_old);
    std::copy(in_edges, in_edges+n_nodes, in_edges_old);
    for(i = 0 ;i<n_nodes;i++){
        int to = reindexedNew[i];
        owner[to] = owner_old[i];
        priority[to] = priority_old[i];
        out_edges[to] = out_edges_old[i];
        in_edges[to] = in_edges_old[i];

        for(auto & j : out_edges[to]){
            j = make_tuple(reindexedNew[target(j)], guard_index(j));
        }
        for(auto & j : in_edges[to]){
            j = make_tuple(reindexedNew[target(j)], guard_index(j));
        }
    }
    for(int & edge_origin : edge_origins){
        edge_origin = reindexedNew[edge_origin];
    }
//    unordered_set<int> p_old;
//    for(auto & p : priorityI){
//        p_old = p;
//        p.clear();
//        for(const auto & pp : p_old)
//            p.insert(reindexedNew[pp]);
//    }
    delete[] out_edges_old;
    delete[] in_edges_old;
}

void Game::parsePGFromFile(const string &filename) {
    parsePG = true;
    parseVPGFromFile(filename);
}

void Game::writePG(ostream *output) {
    *output << "parity " << n_nodes << ';';
    for(int v = 0;v < n_nodes;v++){
        *output << endl << reindexedOrg[v] << ' ' << priority[v] << ' ' << owner[v];
        char seperator = ' ';
        for(const auto & e : out_edges[v]) {
            *output << seperator << reindexedOrg[target(e)];
            seperator = ',';
        }
        *output << ';';
    }
}

void Game::compressVertices() {
    compressvertices = true;
    orgvertices.resize(n_nodes);
    for(int i = 0;i < n_nodes;i++){
        orgvertices[i].resize(1);
        orgvertices[i][0] = i;
    }
    vector<bool> remove(n_nodes);
    for(int i = 0;i < n_nodes;i++){
        int t = target(out_edges[i][0]);
        if(out_edges[i].size() ==  1 && t != i && priority[t] >= priority[i]){
            int orgindex = orgvertices[t].size();
            orgvertices[t].resize(orgindex + orgvertices[i].size());
            for(int j = 0;j < orgvertices[i].size();j++)
                orgvertices[t][orgindex + j] = orgvertices[i][j];
            // assume we admit all configurations since VPGs are total
            remove[i] = true;
            for(auto &e : in_edges[i]){
                bool merge = false;
                tuple<int,int> edgeorg, edgenew;

                for(auto & oe: out_edges[target(e)]){
                    if(target(oe) == t) {
                        merge = true;
                        edgeorg = oe;
                    }
                    if(target(oe) == i){
                        oe = make_tuple(t, guard_index(oe));
                        edgenew = oe;
                    }
                }
                if(merge) {
                    cout << "Merge"<<endl;
                    //todo remove this loop by storing the iterator in the previous loop
                    auto ite = out_edges[target(e)].begin();
                    while(ite != out_edges[target(e)].end()){
                        if(*ite == edgenew){
                            out_edges[target(e)].erase(ite);
                            break;
                        }
                        ite++;
                    }
                    edge_guards[guard_index(edgeorg)] |= edge_guards[guard_index(edgenew)];
                }else{
                    int index = in_edges[t].size();
                    in_edges[t].resize(index + 1);
                    in_edges[t][index] = e;
                }
            }
            auto ite = in_edges[t].begin();
            while(ite != in_edges[t].end()){
                if(target(*ite) == i){
                    in_edges[t].erase(ite);
                    break;
                }
                ite++;
            }
            cout << "Remove " << i << endl;
        }
    }

//    vector<bool> remove(n_nodes);
    vector<int> reindex(n_nodes);
    for(int i = 0;i < n_nodes;i++)
        reindex[i] = i;
    int removed = 0;
    for(int i = 0;i < n_nodes;i++){
        if(remove[i])
            removed++;
        else
            reindex[i] -= removed;
    }
    for(int i = 0;i < n_nodes;i++){
        if(!remove[i]) {
            if(reindex[i] < i){
                owner[reindex[i]] = owner[i];
                priority[reindex[i]] = priority[i];
                out_edges[reindex[i]] = out_edges[i];
                in_edges[reindex[i]] = in_edges[i];
                orgvertices[reindex[i]] = orgvertices[i];
            } else if(reindex[i] > i)
                throw "Incorrect reindex";
        }
    }
    n_nodes -= removed;
    priority.resize(n_nodes);
    owner.resize(n_nodes);
    orgvertices.resize(n_nodes);
    for(int i = 0;i < n_nodes;i++){
        for(auto & e : out_edges[i])
            e = make_tuple(reindex[target(e)], guard_index(e));
        for(auto & e : in_edges[i])
            e = make_tuple(reindex[target(e)], guard_index(e));
    }

//    for(int i = 0;i < n_nodes;i++){
//        int t = reindex[target(out_edges[i][0])];
//        if(out_edges[i].size() ==  1 && t != i){
//            // assume we admit all configurations since VPGs are total
//            remove[i] = true;
//            reindex[i] = t;
//            if(priority[i] > priority[t]){
//                priorityI[priority[t]].erase(t);
//                priorityI[priority[i]].insert(t);
//                priority[t] = priority[i];
//                for(int j = 0;j < n_nodes;j++)
//                    if(reindex[j] == i)
//                        reindex[j] = t;
//            }
//            //todo remove edges
////            out_edges[t].reserve(out_edges[t].size() + out_edges[i].size());
////            out_edges[t].insert(out_edges[t].end(), out_edges[i].begin(), out_edges[i].end());
//        }
//    }

//    auto * out_edges_old = out_edges;
//    out_edges = new vector<std::tuple<int,int>>[n_nodes - removed];
//    vector<int> priority_old(n_nodes);
//    priority_old = priority;
//    vector<int> owner_old(n_nodes);
//    owner_old = owner;
//    for(int i = 0;i < n_nodes;i++){
//        if(!remove[i] && reindex[i] != i){
//            out_edges[reindex[i]] = out_edges_old[i];
//            priority[reindex[i]] = priority_old[i];
//            owner[reindex[i]] = owner_old[i];
//        }
//    }
//    delete[] out_edges_old;
//
//    orgvertices.resize(n_nodes - removed);
//    for(int i = 0;i < n_nodes;i++){
//        orgvertices[reindex[i]].insert(i);
//    }
//    n_nodes -= removed;
//    //todo vector instead of array of out_edges
////    delete[] (out_edges + n_nodes);
//
//
//    for(auto &p : priorityI){
//        unordered_set<int> pp = p;
//        p.clear();
//        for(const auto &v : pp)
//            if (!remove[v])
//                p.insert(reindex[v]);
//    }
//    //todo remove this reset
//    for(int i = 0;i < n_nodes;i++)
//        in_edges[i].resize(0);
//    buildInEdges();
}

void Game::buildInEdges() {
    for(int i = 0;i < n_nodes;i++){
        for(const auto & e : out_edges[i]){
            int t = target(e);
            int index = in_edges[t].size();
            in_edges[t].resize(index + 1);
            in_edges[t][index] = make_tuple(i, guard_index(e));
        }
    }
}
