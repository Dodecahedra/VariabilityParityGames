#include <iostream>
#include <chrono>
#include "Game.h"
#include "zlnk.h"
#include "FPIte.h"
#include "MBR.h"

using namespace std;

string winningset;
void printFull(char* varset, int size, string prod){
    if(size == 0){
        cout << "For product " << prod << " the following vertices are in: " << winningset << endl;
        return;
    }

    if(*varset < 0){
        string prod2 = prod;
        prod.append("0");
        prod2.append("1");
        printFull(varset + 1, size -1, prod);
        printFull(varset + 1, size -1, prod2);
    } else {
        if(*varset == 0){
            prod.append("0");
        } else {
            prod.append("1");
        }
        printFull(varset +1, size -1, prod);
    }
}
void allsatPrintHandler(char* varset, int size)
{
    string p = "";
    printFull(varset, size, p);
}

int main(int argc, char** argv) {
#ifdef SINGLEMODE
    cout << "== Running Single mode == \n";
#endif
#ifdef subsetbdd
    cout << "== Subsets are bdds == \n";
#endif
#ifdef subsetexplicit
    cout << "== Subsets are explicit == \n";
#endif
    if(argc < 2) {
        cerr << "Incorect params";
        return 2;
    }
    try {
        time_t t0 = time(0);
        cout << "\n[" << t0 << "] Start\n";

        Game g;
        char *specificconf;
        bool specificconfenabled = false;

        bool fulloutput = false;
        bool metricoutput = false;
        bool priocompress = false;
        bool SolveFPIte = false;

        for(int i = 2;i<argc;i++){
            switch (*argv[i]){
                case 'c':
                    specificconfenabled = true;
                    specificconf = argv[i]+1;
                    break;
                case 'f':
                    fulloutput = true;
                    break;
                case 'm':
                    metricoutput = true;
                    break;
                case 'p':
                    priocompress = true;
                    break;
                case 'F':
                    SolveFPIte = true;
                    priocompress = true;
                    break;
                default:
                    cerr << "Unknown parameter: " << argv[i];
                    break;
            }
        }
        if(specificconfenabled){
            g.parseGameFromFile(argv[1],specificconf);
        } else {
            g.parseGameFromFile(argv[1]);
        }
        time_t t = time(0);

        cout << "\n[" << t << "] Parsed\n";
        fflush(stdout);

#ifdef subsetbdd
        // enable cache after parsing
        bdd_setcacheratio(200);
#endif

        auto start = std::chrono::high_resolution_clock::now();

        if(priocompress)
            g.compressPriorities();

        if(SolveFPIte){
#ifdef SINGLEMODE
            FPIte fpite(&g);
            fpite.solve();
            auto end = std::chrono::system_clock::now();

            auto elapsed =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            cout << "Solving time: " << elapsed.count() << " ns";
            time_t t2 = time(0);
            cout << '[' << t2 << "] Solved\n";

            cout << "W0: \n";
            cout << "The following vertices are in: ";
            for(int i = 0;i<(*fpite.W0).size();i++){
                if((*fpite.W0)[i])
                    cout << i << ",";
            }
            cout << "\nW1: \n";
            cout << "The following vertices are in: ";
#else
            MBR mbr(&g);
            mbr.solve();
            auto end = std::chrono::system_clock::now();

            auto elapsed =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            cout << "Solving time: " << elapsed.count() << " ns";
            time_t t2 = time(0);
            cout << '[' << t2 << "] Solved\n";

#ifdef subsetbdd
            cout << "W0: \n";
            for(int i =0;i<MBR::winningConf.size();i++) {
                winningset = "";
                if(fulloutput){
                    for (int j = 0; j < MBR::winningVertices[i].size(); j++){
                        if (MBR::winningVertices[i][j]) {
                            winningset.append(std::to_string(j));
                            winningset.append(",");
                        }
                    }
                } else if(MBR::winningVertices[i][0]){
                    winningset = "0,";
                }
                bdd_allsat(MBR::winningConf[i], allsatPrintHandler);
            }
            cout << "W1: \n";
            for(int i =0;i<MBR::winningConf.size();i++) {
                winningset = "";
                if(fulloutput){
                    for (int j = 0; j < MBR::winningVertices[i].size(); j++){
                        if (!MBR::winningVertices[i][j]) {
                            winningset.append(std::to_string(j));
                            winningset.append(",");
                        }
                    }
                } else if(!MBR::winningVertices[i][0]){
                    winningset = "0,";
                }
                bdd_allsat(MBR::winningConf[i], allsatPrintHandler);
            }
#else
            cout << "Can't print output";
#endif
#endif
        } else {
            zlnk::conf_metricoutput = metricoutput;
            zlnk z(&g);
            auto * W0BigV = new unordered_set<int>;
            auto * W1BigV = new unordered_set<int>;

    #ifdef SINGLEMODE
            vector<Subset> * W0vc = nullptr;
            vector<Subset> * W1vc = nullptr;
    #else
            vector<Subset> * W0vc = new vector<Subset>(g.n_nodes);
            vector<Subset> * W1vc = new vector<Subset>(g.n_nodes);
            for(int i = 0;i<g.n_nodes;i++){
                (*W0vc)[i] = emptyset;
                (*W1vc)[i] = emptyset;
            }
    #endif

            z.solve(W0BigV, W0vc, W1BigV, W1vc);


            auto end = std::chrono::system_clock::now();

            auto elapsed =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            cout << "Solving time: " << elapsed.count() << " ns";
            time_t t2 = time(0);
            cout << '[' << t2 << "] Solved\n";

            cout << "W0: \n";
            g.printCV(W0BigV, W0vc, fulloutput);
            cout << "W1: \n";
            g.printCV(W1BigV, W1vc, fulloutput);
        }
    } catch(std::string s)
    {
        cerr << s << "\n";
        return 1;
    }
    return 0;
}
