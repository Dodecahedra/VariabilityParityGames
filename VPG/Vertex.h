//
// Created by sjef on 5-6-19.
//

#ifndef VPGSOLVER_VERTEX_H
#define VPGSOLVER_VERTEX_H


#include "Edge.h"
#include <vector>

class Edge;
class Vertex {
public:
    int owner;
    int priority;
    Vertex();
    Vertex(int owner, int priority);
    std::vector<Edge*> out;
    std::vector<Edge*> in;
    void addOutgoingEdge(Edge *e);
    void removeOutgoingEdge(Edge *e);
    void removeOutgoingEdge(int i);
    void removeOutgoingEdge(std::vector<Edge*>::iterator iter, Edge *e);


protected:
    void removeIncomingEdge(Edge *e);
};


#endif //VPGSOLVER_VERTEX_H
