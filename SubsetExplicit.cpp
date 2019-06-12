//
// Created by sjef on 11-6-19.
//

#include <algorithm>
#include "SubsetExplicit.h"

SubsetExplicit SubsetExplicit::SetFullset = SubsetExplicit();
SubsetExplicit SubsetExplicit::SetEmptyset = SubsetExplicit();

SubsetExplicit *SubsetExplicit::getFullSet() {
    return &SubsetExplicit::SetEmptyset;
}

SubsetExplicit *SubsetExplicit::getEmptySet() {
    return &SubsetExplicit::SetEmptyset;
}


void SubsetExplicit::operator=(SubsetExplicit other) {
    this->items =  other.items;
}

void SubsetExplicit::operator&=(SubsetExplicit &other) {
    auto thisite = this->items.begin();
    auto otherite = other.items.begin();

    while(thisite != this->items.end()){
        if(otherite == other.items.end() || *thisite < *otherite) {
            thisite = this->items.erase(thisite);
        } else if(*thisite == *otherite) {
            thisite++;
            otherite++;
        } else {
            //*thisite > *otherite
            otherite++;
        }
    }
}

void SubsetExplicit::operator|=(SubsetExplicit &other) {
    // todo improve
    for(const auto &i : other.items){
        this->items.insert(i);
    }
}

void SubsetExplicit::operator-=(SubsetExplicit &other) {
    auto thisite = this->items.begin();
    auto otherite = other.items.begin();

    while(thisite != this->items.end() && otherite != other.items.end()){
        if(*thisite == *otherite) {
            thisite = this->items.erase(thisite);
            otherite++;
        } else if(*thisite < *otherite){
            thisite++;
        } else {
            //*thisite > *otherite
            otherite++;
        }
    }
}

SubsetExplicit::SubsetExplicit() {

}

int SubsetExplicit::operator==(SubsetExplicit &other) {
    return this->items == other.items;
}

SubsetExplicit::SubsetExplicit(int range, int bit) {
    int mask = 1 << bit;
    for(int i = 0;i<=range;i++)
        if((i & mask) != 0)
            this->items.insert(i);
}
