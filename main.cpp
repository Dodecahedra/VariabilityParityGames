#include <iostream>
#include <chrono>
#include "Game.h"
#include "zlnk.h"

using namespace std;

int main(int argc, char** argv) {
    if(argc < 2) {
        cerr << "Incorect params";
        return 2;
    }
    try {
        time_t t0 = time(0);
        cout << "\n[" << t0 << "] Start\n";

        Game g;
        g.parseGameFromFile(argv[1]);
        time_t t = time(0);
        cout << "\n[" << t << "] Parsed\n";
        fflush(stdout);

        auto start = std::chrono::high_resolution_clock::now();


        zlnk z(&g);
        auto * W0BigV = new unordered_set<int>;
        vector<bdd> * W0vc = new vector<bdd>(g.n_nodes);

        auto * W1BigV = new unordered_set<int>;
        vector<bdd> * W1vc = new vector<bdd>(g.n_nodes);

        for(int i = 0;i<g.n_nodes;i++){
            (*W0vc)[i] = bddfalse;
            (*W1vc)[i] = bddfalse;
        }
        z.solve(W0BigV, W0vc, W1BigV, W1vc);

        auto end = std::chrono::system_clock::now();

        auto elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        cout << "Solving took " << elapsed.count() << "ms.\n";
        time_t t2 = time(0);
        cout << '[' << t2 << "] Solved\n";
        cout << "W0: \n";
        g.printCV(W0BigV, W0vc);
        cout << "W1: \n";
        g.printCV(W1BigV, W1vc);
    } catch(std::string s)
    {
        cerr << s << "\n";
        return 1;
    }
    return 0;
}