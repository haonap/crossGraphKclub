//
// Created by Hao Pan on 9/20/21.
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
    vector<int> adjacency; // used in CommunityPeel
    vector<int> F; // used in CommunityPeel
    vector<int> numCommonNeighbors; // used in CommunityPeel
    vector<vector<int>> commonNeighbors; // used in CommunityPeel

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
    void initialize();// initialize graph, set #vertex = n, and all others 0 or empty
    vector<int> FindCommonNeighbors(int, int);
    vector<int>* FindCommonKNeighbors(int, int);
    vector<int> KBFS(int, int);
    void DeleteEdge(int u, int v); // remove edge uv from graph: make sure uv is an edge
    graph* CreateSubgraph(vector<int>* C);//create a subgraph induced by set C
    //graph remove nodes C: if node u in removedC, then removedC[u] = true; the function will reset the node degree =0 and removed it from its neighbors
    void GraphRemoveNodes(const vector<bool> &removedC);
    vector<int> BFS(int);
    void CalculateDistanceFrom(int);
    void CalculateDistanceFromTo(int, vector<bool>*);
    vector<int> FindDegeneracyOrdering(vector<int>&); // from parsimonious paper
    vector<int> FindHeuristicClique(vector<int>&, vector<int>&); // from parsimonious paper
    vector<int> ShortestPathsUnweighted(int, vector<bool>&); // from parsimonious paper
    vector<int> MultiSourceShortestPaths(vector<int>&, vector<bool>&); // from parsimonious paper
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
    vector<GRBModel>& SEP_MODELS;
    vector<GRBVar*>& Z;
    vector<GRBVar*>& W;
    double epsilon;
    int& independentSetInequalityCount;
    int& independentSetInequalityCountViolateEpsilon;

    pGraph_callback(string, GRBVar*, vector<graph*>&, int, int&, int&, int&, double&, double&, double&, vector<GRBModel>&, vector<GRBVar*>&, vector<GRBVar*>&, double, int&, int&);

protected:
    //callback function
    void callback();
};


//callback for Independent set equality to check if node is at the root node: for experiment test
class mycallback : public GRBCallback{
public:
    double lastiter;
    double lastnode;
    int numvars;
    GRBVar* vars;
    mycallback(GRBVar* xvars) {
        vars = xvars;
    }

    //mycallback();
protected:
    //callback function
    void callback(){
        try {
            if (where == GRB_CB_POLLING) {
                //cout<<"i am in callback"<<endl;
                // Ignore polling callback
            } else if (where == GRB_CB_MIPSOL) {
                //cout<<"i am in callback MIP SOL"<<endl;
                // MIP solution callback
                int nodecnt = (int) getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
                int solcnt = getIntInfo(GRB_CB_MIPSOL_SOLCNT);
                if (nodecnt >= 1) {
                    cout << "the number of node explored is: " << nodecnt << endl;
                    cout << "the solution status is:" << solcnt << endl;
                    cout << "Stop early - root node explored" << endl;
                    abort();
                }
            }
        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};

#endif //PGRAPHKCLUB_CLASSES_H



