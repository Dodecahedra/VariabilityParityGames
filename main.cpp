/***********************************************************************************************************************
 * In this project only algorithms and datastructures are implemented properly.
 *
 * The code in this file is not optimized and not up to standards; it is sufficient only for experimental applications
 * but not for any real application.
 **********************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <chrono>
#include "VPGame.h"
#include "Algorithms/FPIte.h"
#include "Algorithms/MBR.h"

#ifdef SINGLEMODE
#include "Algorithms/zlnkPG.h"
#else
#include "Algorithms/zlnkVPG.h"
#endif

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
    if (argc < 2) {
        cerr << "Incorect params";
        return 2;
    }
    try {
        time_t t0 = time(0);
        cout << "\n[" << t0 << "] Start\n";

        VPGame g;

        char *specificconf;
        bool specificconfenabled = false;

        bool fulloutput = false;
        bool metricoutput = false;
        bool priocompress = false;
        bool SolveFPIte = false;
        bool assistedW0 = false;
        char *assistanceW0;
        bool assistrandom = false;
        int assistn;
        bool solvelocal = false;
        bool regulargame = false;
        bool compressvertices = false;
        char *metricdir;
        bool projectmode = false;
        char *projectdir;
        for (int i = 2; i < argc; i++) {
            switch (*argv[i]) {
                case 'Q':
                    projectmode = true;
                    projectdir = argv[i] + 1;
                    break;
                case 'v':
                    compressvertices = true;
                    break;
                case 'o':
                    g.specificvarlast = true;
                    g.specificvar = atoi(argv[i] + 1);
                    break;
                case 'R':
                    regulargame = true;
                    break;
                case 'l':
                    solvelocal = true;
                    break;
                case 'P':
                    assistedW0 = true;
                    assistanceW0 = argv[i] + 1;
                    break;
                case 'a':
                    assistrandom = true;
                    assistn = atoi(argv[i] + 1);
                    break;
                case 'c':
                    specificconfenabled = true;
                    specificconf = argv[i] + 1;
                    break;
                case 'f':
                    fulloutput = true;
                    break;
                case 'm':
                    metricoutput = true;
                    metricdir = argv[i] + 1;
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
        vector<int> loops = vector<int>();
        if (regulargame) {
            g.parsePGFromFile(argv[1]);
        } else {
            if (specificconfenabled) {
                g.parseVPGFromFile(argv[1], specificconf);
            } else {
                g.parseVPGFromFile(argv[1]);
            }
        }
        if (projectmode) {
            vector<tuple<ConfSet, string>> allconfs;
            g.findAllElements(g.bigC, &allconfs);
            for (const auto &conf : allconfs) {
                string target = string(projectdir);
                target += "/sSVPG" + get<1>(conf);
                ofstream output;
                output.open(target);
                g.writePG(&output, get<0>(conf));
                output.close();
                cout << "Projected to " << target << endl;
            }
            return 0;
        }
#ifdef lambdameasure
        float avg = 0;
        int c = 0;
        for(const auto & e : g.edge_guards){
            float n =  bdd_satcount(e);
            cout << "Edge allows " << n << " configurations." << endl;
            avg += (n - avg) / ++c;
        }
        cout << "Edge allow on average " << avg << " configurations." << endl;
        return 0;
#endif
        time_t t = time(0);

        cout << "\n[" << t << "] Parsed\n";
        fflush(stdout);

#ifdef subsetbdd
        bdd_gbc();
        // enable cache after parsing
        bdd_setcacheratio(200);
        bdd_gbc();
#endif
        auto start = std::chrono::high_resolution_clock::now();

        if(compressvertices)
            g.compressVertices();

        if(priocompress)
            g.compressPriorities();

        if(SolveFPIte){
            g.reindexVertices();
        start = std::chrono::high_resolution_clock::now();
#ifdef SINGLEMODE
            FPIte fpite(&g);
            fpite.solvelocal = solvelocal;
            if(assistedW0)
                fpite.setP0(assistanceW0);
            if(assistrandom){
                fpite.P0IsFull();
                fpite.unassist(g.n_nodes - assistn);
            }
            fpite.solve();
            auto end = std::chrono::high_resolution_clock::now();

            auto elapsed =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            cout << "Solving time: " << elapsed.count() << " ns";
            time_t t2 = time(0);
            cout << '[' << t2 << "] Solved\n";
            cout << "Executed db " << fpite.dbs_executed << ',' << fpite.verticesconsidered << " vertices" << endl;
            cout << "W0: \n";
            cout << "The following vertices are in: ";
            if(fulloutput)
            {
                for(int i = 0;i<(*fpite.W0).size();i++){
                    if((*fpite.W0)[i]) {
                        for(const auto & v : g.orgvertices[g.reindexedOrg[i]]) {
                            cout << v << ",";
                        }
                    }
                }
            } else {
                if((*fpite.W0)[g.reindexedNew[g.findVertexWinningFor0()]])
                    cout << "0,";
            }
            cout << "\nW1: \n";
            cout << "The following vertices are in: ";
#else
            MBR mbr(&g);
            mbr.solvelocal = solvelocal;
            MBR::fulloutput = fulloutput;
            MBR::metric_output = metricoutput;
            if(metricoutput)
                MBR::metric_dir = metricdir;
            mbr.solve();
            auto end = std::chrono::high_resolution_clock::now();

            auto elapsed =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            cout << "Solving time: " << elapsed.count() << " ns";
            time_t t2 = time(0);
            cout << '[' << t2 << "] Solved\n";
            if(metricoutput)
                mbr.printMeasurements(&cout);
#ifdef subsetbdd
            cout << "W0: \n";
            for(int i =0;i<MBR::winningConf.size();i++) {
                winningset = "";
                if(fulloutput){
                    for (int j = 0; j < MBR::winningVertices[i].size(); j++){
                        if (MBR::winningVertices[i][j]) {
                            for(const auto & v : g.orgvertices[g.reindexedOrg[j]]) {
                                winningset.append(std::to_string(v));
                                winningset.append(",");
                            }
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
                            for(const auto & v : g.orgvertices[g.reindexedOrg[j]]) {
                                winningset.append(std::to_string(v));
                                winningset.append(",");
                            }
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
            auto * W0BigV = new VertexSetZlnk(g.n_nodes);
            auto * W1BigV = new VertexSetZlnk(g.n_nodes);

    #ifdef SINGLEMODE
            zlnkPG::conf_metricoutput = metricoutput;
            zlnkPG z(&g);

            vector<ConfSet> * W0vc = nullptr;
            vector<ConfSet> * W1vc = nullptr;
            if(solvelocal)
                z.solvelocal = 2;
            z.solve(W0BigV, W1BigV);
    #else
            zlnkVPG::conf_metricoutput = metricoutput;
            zlnkVPG z(&g);
            vector<ConfSet> * W0vc = new vector<ConfSet>(g.n_nodes);
            vector<ConfSet> * W1vc = new vector<ConfSet>(g.n_nodes);
            for(int i = 0;i<g.n_nodes;i++){
                (*W0vc)[i] = emptyset;
                (*W1vc)[i] = emptyset;
            }
            if(solvelocal)
                z.solvelocal = 2;
            z.solve(W0BigV, W0vc, W1BigV, W1vc);
    #endif


            auto end = std::chrono::high_resolution_clock::now();

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