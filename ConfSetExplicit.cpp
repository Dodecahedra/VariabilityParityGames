//
// Created by sjef on 11-6-19.
//

#include <algorithm>
#include "ConfSetExplicit.h"


ConfSetExplicit ConfSetExplicit::SetFullset = ConfSetExplicit();
ConfSetExplicit ConfSetExplicit::SetEmptyset = ConfSetExplicit();
int ConfSetExplicit::size;


ConfSetExplicit& ConfSetExplicit::operator=(ConfSetExplicit other) {
    this->items = other.items;
    return *this;
}

void ConfSetExplicit::operator&=(ConfSetExplicit &other) {
    this->items &= other.items;
}

void ConfSetExplicit::operator|=(ConfSetExplicit &other) {
    this->items |= other.items;
}

void ConfSetExplicit::operator-=(ConfSetExplicit &other) {
    this->items -= other.items;
}

ConfSetExplicit::ConfSetExplicit() {
    items.resize(size);
}

int ConfSetExplicit::operator==(ConfSetExplicit &other) {
    return this->items == other.items;
}

ConfSetExplicit::ConfSetExplicit(int bit) {
    items.resize(size);
    unsigned int mask = 1 << bit;
    for(unsigned int i = 0;i<=ConfSetExplicit::size;i++)
        this->items[i] = (i & mask) != 0;
}

int ConfSetExplicit::count() {
    int c = 0;
    for(auto b : items)
        if(b)
            c++;
    return c;
}
