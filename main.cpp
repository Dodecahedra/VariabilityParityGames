#include <iostream>
#include "Game.h"

using namespace std;

int main(int argc, char** argv) {
    if(argc < 2) {
        cerr << "Incorect params";
        return 2;
    }
    try {

        Game g;
        g.parseGameFromFile(argv[1]);

    } catch(std::string s)
    {
        cerr << s << "\n";
        return 1;
    }
    return 0;
}