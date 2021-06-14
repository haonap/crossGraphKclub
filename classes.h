//
// Created by Hao Pan on 4/5/21.
//

#ifndef PGRAPHKCLUB_CLASSES_H
#define PGRAPHKCLUB_CLASSES_H


#include "gurobi_c++.h"
#include "functions.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
using namespace std;


class vertex{
public:
    int ind;
    int name;
    int degree;
    vector<int> neighbors;
    vector<int> kNeighbors;
    vector<int> distanceTo;

    vertex(int);
    vertex(){};
};


class graph {
public:
    int ind;
    int n;
    int m;
    int maxDegree;
    int maxDegreeVertexIndex;
    vector<vertex*> vertices;

    graph();
    ~graph();

    bool IsAdjacent(int, int);
    bool IsKAdjacent(int, int);
    void FindKNeighbors(int);
    vector<int>* FindCommonKNeighbors(int, int);
    vector<int> KBFS(int, int);
    graph* CreateSubgraph(vector<int>*);
    void CalculateDistanceFrom(int);
    void CalculateDistanceFromTo(int, vector<bool>*);
};


class pGraph_callback : public GRBCallback{
public:

    GRBVar* X;
    vector<graph*>& graphCollection;
    string methodParam;
    int kParam;
    int& numStrengthenedLazyCut;
    int& numLazyCut;
    int& numCallback;
    double& callbackTime;
    double& aveNumTermsNeg;
    double& aveNumTermsNegStrengthened;

    pGraph_callback(string, GRBVar*, vector<graph*>&, int, int&, int&, int&, double&, double&, double&);

protected:
    //callback function
    void callback();
};


#endif //PGRAPHKCLUB_CLASSES_H
