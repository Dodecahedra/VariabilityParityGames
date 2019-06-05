//
// Created by sjef on 5-6-19.
//

#include <algorithm>
#include <vector>
#include "Vertex.h"

Vertex::Vertex() {
    Vertex(0,0);
}

Vertex::Vertex(int owner, int priority) {
    this->owner = owner;
    this->priority = priority;
}

void Vertex::addOutgoingEdge(Edge *e) {
    int n = out.size();
    out.resize(n + 1);
    out[n] = e;
    n = e->target->in.size();

    e->target->in.resize(n + 1);
    e->target->in[n] = e;
}

void Vertex::removeOutgoingEdge(Edge *e) {
    auto iter = std::find(out.begin(), out.end(), e);
    if(iter != out.end())
    {
        removeOutgoingEdge(iter, e);
    }
}

void Vertex::removeOutgoingEdge(std::vector<Edge *>::iterator iter, Edge* e) {
    out.erase(iter);
    e->target->removeIncomingEdge(e);
    delete e;
}

void Vertex::removeOutgoingEdge(int i) {
    removeOutgoingEdge(out.begin() + i, out[i]);
}

void Vertex::removeIncomingEdge(Edge *e) {
    auto iter = std::find(in.begin(), in.end(), e);
    if(iter != in.end())
    {
        in.erase(iter);
    }
}
