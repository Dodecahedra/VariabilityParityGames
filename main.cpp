#include <iostream>
#include "Game.h"
#include "zlnk.h"

using namespace std;

int main(int argc, char** argv) {
    if(argc < 2) {
        cerr << "Incorect params";
        return 2;
    }
    try {

        Game g;
        g.parseGameFromFile(argv[1]);
        zlnk z(&g);
        auto * W0BigV = new unordered_set<int>;
        vector<BDD> * W0vc = new vector<BDD>(g.n_nodes);

        auto * W1BigV = new unordered_set<int>;
        vector<BDD> * W1vc = new vector<BDD>(g.n_nodes);

        for(int i = 0;i<g.n_nodes;i++){
            (*W0vc)[i] = g.bm->getZero();
            (*W1vc)[i] = g.bm->getZero();
        }
        time_t t = time(0);
        cout << "\n[" << t << "] Parsed\n";
        fflush(stdout);
        z.solve(W0BigV, W0vc, W1BigV, W1vc);
        time_t t2 = time(0);
        cout << '[' << t << "] Solved\n";
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