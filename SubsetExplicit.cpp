//
// Created by sjef on 11-6-19.
//

#include <algorithm>
#include "SubsetExplicit.h"


SubsetExplicit SubsetExplicit::SetFullset = SubsetExplicit();
SubsetExplicit SubsetExplicit::SetEmptyset = SubsetExplicit();
int SubsetExplicit::size;

SubsetExplicit *SubsetExplicit::getFullSet() {
    return &SubsetExplicit::SetFullset;
}

SubsetExplicit *SubsetExplicit::getEmptySet() {
    return &SubsetExplicit::SetEmptyset;
}


void SubsetExplicit::operator=(SubsetExplicit other) {
    this->items = other.items;
}

void SubsetExplicit::operator&=(SubsetExplicit &other) {
    this->items &= other.items;
}

void SubsetExplicit::operator|=(SubsetExplicit &other) {
    this->items |= other.items;
}

void SubsetExplicit::operator-=(SubsetExplicit &other) {
    this->items -= other.items;
}

SubsetExplicit::SubsetExplicit() {
    items.resize(size);
}

int SubsetExplicit::operator==(SubsetExplicit &other) {
    return this->items == other.items;
}

SubsetExplicit::SubsetExplicit(int bit) {
    items.resize(size);
    int mask = 1 << bit;
    for(int i = 0;i<=SubsetExplicit::size;i++)
        this->items[i] = (i & mask) != 0;
}

int SubsetExplicit::count() {
    int c = 0;
    for(auto b : items)
        if(b)
            c++;
    return c;
}
