//
// Created by Hao Pan on 4/5/21.
//


#include <iostream>
#include "functions.h"
#include "classes.h"
using namespace std;

graph::graph(){};

graph::~graph(){

    for(auto vertexPtr : vertices)
        delete vertexPtr;
}


bool graph::IsAdjacent(int u, int v) {

    if(binary_search(vertices[u]->neighbors.begin(), vertices[u]->neighbors.end(), v))
        return true;
    return false;
}


bool graph::IsKAdjacent(int u, int v) {

    if(binary_search(vertices[u]->kNeighbors.begin(), vertices[u]->kNeighbors.end(), v))
        return true;
    return false;
}


void graph::FindKNeighbors(int k){

    for(int i = 0; i < n; i++){
        vertices[i]->kNeighbors.clear();
        vertices[i]->kNeighbors = KBFS(i, k);
    }
}

vector<int>* graph::FindCommonKNeighbors(int u, int v){

    vector<int>* commonKNeighbors = new vector<int>();
    commonKNeighbors->clear();
    int i = 0, j = 0;
    while(i < vertices[u]->kNeighbors.size() && j < vertices[v]->kNeighbors.size()){
        if(vertices[u]->kNeighbors[i] == vertices[v]->kNeighbors[j]){
            commonKNeighbors->push_back(vertices[u]->kNeighbors[i]);
            i++;
            j++;
        }else if(vertices[u]->kNeighbors[i] > vertices[v]->kNeighbors[j]){
            j++;
        }else{
            i++;
        }

    }
    return commonKNeighbors;
}


vector<int> graph::KBFS(int s, int k){

    vector<int> tempVec;
    vector<int>* distance = new vector<int>(n, n);
    queue<int>* Q = new queue<int>;
    (*distance)[s] = 0;
    Q->push(s);
    int u, v, goFlag = 1;

    while(Q->size() && goFlag){
        u = Q->front();
        Q->pop();
        for(int i = 0; i < vertices[u]->neighbors.size(); i++){
            v = vertices[u]->neighbors[i];
            if((*distance)[v] > (*distance)[u] + 1){
                (*distance)[v] = (*distance)[u] + 1;
                if((*distance)[v] > k){
                    goFlag = 0;
                    break;
                }else{
                    Q->push(v);
                }
            }
        }

    }

    for(int i = 0; i < n; i++){
        if((*distance)[i] <= k && i != s)
            tempVec.push_back(i);
    }

    delete distance;
    delete Q;

    return tempVec;
}


graph* graph::CreateSubgraph(vector<int>* C){

    graph* graphPtr = new graph();
    graphPtr->n = C->size();

    for(int i = 0; i < C->size(); i++){
        vertex* vertexPtr = new vertex();
        vertexPtr->ind = i;
        vertexPtr->name = (*C)[i];
        graphPtr->vertices.push_back(vertexPtr);
    }

    graphPtr->m = 0;
    for(int i = 0; i < graphPtr->n; i++){
        int p = 0, q = 0, pMax = vertices[(*C)[i]]->degree, qMax = C->size();
        while(p < pMax && q < qMax){
            int u = vertices[(*C)[i]]->neighbors[p], v = (*C)[q];
            if(u == v){
                graphPtr->vertices[i]->neighbors.push_back(q);
                graphPtr->vertices[q]->name = v;
                p++, q++;
            }else if(u > v)
                q++;
            else
                p++;
        }
        graphPtr->vertices[i]->degree = graphPtr->vertices[i]->neighbors.size();
        graphPtr->m += graphPtr->vertices[i]->degree;
    }

    graphPtr->m /= 2;
    return graphPtr;
}

void graph::CalculateDistanceFrom(int source){

    vector<bool>* S = new vector<bool>(n, true);
    CalculateDistanceFromTo(source, S);
}

void graph::CalculateDistanceFromTo(int source, vector<bool>* S) {

    queue<int>* Q = new queue<int>();
    vertices[source]->distanceTo.clear();
    vertices[source]->distanceTo.resize(n, n);
    vertices[source]->distanceTo[source] = 0;
    Q->push(source);
    int u, v;

    while(Q->size()){
        u = Q->front();
        Q->pop();
        for(int i = 0; i < vertices[u]->degree; i++){
            v = vertices[u]->neighbors[i];
            if((*S)[v] == false) continue;
            if(vertices[source]->distanceTo[v] > vertices[source]->distanceTo[u] + 1){
                vertices[source]->distanceTo[v] = vertices[source]->distanceTo[u] + 1;
                Q->push(v);
            }
        }
    }

    delete Q;
}


pGraph_callback::pGraph_callback(string methodParam, GRBVar* X, vector<graph*>& graphCollection, int kParam, int& numStrengthenedLazyCut, int& numLazyCut, int& numCallback, double& callbackTime, double& aveNumTermsNeg, double& aveNumTermsNegStrengthened) : methodParam(methodParam), X(X), graphCollection(graphCollection), kParam(kParam), numStrengthenedLazyCut(numStrengthenedLazyCut), numLazyCut(numLazyCut), numCallback(numCallback), callbackTime(callbackTime), aveNumTermsNeg(aveNumTermsNeg), aveNumTermsNegStrengthened(aveNumTermsNegStrengthened){}

void pGraph_callback::callback(){
    try{

        if(where != GRB_CB_MIPSOL) return;

        numCallback++;

        double timeBegin = GetWallTime();

        int nParam = graphCollection[0]->n;
        double* x = getSolution(X, nParam);
        vector<int>* curSol = new vector<int>();
        vector<int>* curSolC = new vector<int>();


        for(int i = 0; i < nParam; i++){
            if(x[i] > 0.5) curSol->push_back(i);
            else curSolC->push_back(i);
        }

        if(curSol->empty()) return;

        graph* subgraphPtr;
        for(auto graphPtr : graphCollection){
            subgraphPtr = graphPtr->CreateSubgraph(curSol);
            for(int i = 0; i < subgraphPtr->n; i++){
                subgraphPtr->CalculateDistanceFrom(i);
                for(int j = i + 1; j < subgraphPtr->n; j++){
                    if(subgraphPtr->vertices[i]->distanceTo[j] > kParam){

                        if(methodParam == "BASE"){

                            vector<int>* minimalSeparatorPtr = MINIMALIZE(graphPtr, subgraphPtr->vertices[i]->name, subgraphPtr->vertices[j]->name, kParam, curSolC);
                            GRBLinExpr EXPR = X[subgraphPtr->vertices[i]->name] + X[subgraphPtr->vertices[j]->name];
                            for(int k = 0; k < minimalSeparatorPtr->size(); k++)
                                EXPR -= X[(*minimalSeparatorPtr)[k]];
                            addLazy(EXPR <= 1);

                            if(numLazyCut == 0) aveNumTermsNeg = minimalSeparatorPtr->size();
                            else aveNumTermsNeg = (aveNumTermsNeg*numLazyCut + minimalSeparatorPtr->size())/(numLazyCut + 1);

                            delete minimalSeparatorPtr;

                        }else if(methodParam == "PGRAPH"){

                            vector<int>* minimalSeparatorPtr = MINIMALIZE(graphPtr, subgraphPtr->vertices[i]->name, subgraphPtr->vertices[j]->name, kParam, curSolC);
                            GRBLinExpr EXPR = X[subgraphPtr->vertices[i]->name] + X[subgraphPtr->vertices[j]->name];
                            bool strengthenedFlag = false;
                            int countNumTermsNeg = 0;
                            for(int k = 0; k < minimalSeparatorPtr->size(); k++){
                                int kVertex = (*minimalSeparatorPtr)[k];
                                bool skipFlag = false;
                                for(auto gPtr : graphCollection){
                                    bool temp1 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[i]->name);
                                    bool temp2 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[j]->name);
                                    if(temp1 && temp2) continue;
                                    else{
                                        skipFlag = true;
                                        break;
                                    }
                                }
                                if(skipFlag == false){
                                    EXPR -= X[kVertex];
                                    countNumTermsNeg++;
                                }
                                else strengthenedFlag = true;
                            }
                            addLazy(EXPR <= 1);

                            if(numLazyCut == 0) aveNumTermsNeg = countNumTermsNeg;
                            else aveNumTermsNeg = (aveNumTermsNeg*numLazyCut + countNumTermsNeg)/(numLazyCut + 1);

                            if(strengthenedFlag == true){
                                if(numStrengthenedLazyCut == 0) aveNumTermsNegStrengthened = countNumTermsNeg;
                                else aveNumTermsNegStrengthened = (aveNumTermsNegStrengthened*numStrengthenedLazyCut + countNumTermsNeg)/(numStrengthenedLazyCut + 1);
                                numStrengthenedLazyCut++;
                            }
                            delete minimalSeparatorPtr;

                        }else if(methodParam == "PGRAPHDELETE"){

                            vector<int>* curSolCPrime = new vector<int>();
                            vector<int>* toBeDeleted = new vector<int>();
                            for(int k = 0; k < curSolC->size(); k++){
                                int kVertex = (*curSolC)[k];
                                bool skipFlag = false;
                                for(auto gPtr : graphCollection){

                                    bool temp1 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[i]->name);
                                    bool temp2 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[j]->name);

                                    if(temp1 && temp2) continue;
                                    else{
                                        skipFlag = true;
                                        break;
                                    }
                                }

                                if(skipFlag == false) curSolCPrime->push_back(kVertex);
                                else toBeDeleted->push_back(kVertex);
                            }

                            GRBLinExpr EXPR = X[subgraphPtr->vertices[i]->name] + X[subgraphPtr->vertices[j]->name];
                            vector<int>* minimalSeparatorPtr = MINIMALIZE(graphPtr, subgraphPtr->vertices[i]->name, subgraphPtr->vertices[j]->name, kParam, toBeDeleted, curSolCPrime);

                            for(auto k : (*minimalSeparatorPtr))
                                EXPR -= X[k];
                            addLazy(EXPR <= 1);

                            if(numLazyCut == 0) aveNumTermsNeg = minimalSeparatorPtr->size();
                            else aveNumTermsNeg = (aveNumTermsNeg*numLazyCut + minimalSeparatorPtr->size())/(numLazyCut + 1);

                            vector<bool>* checkStrength = new vector<bool>(graphPtr->n, true);
                            for(auto k : (*minimalSeparatorPtr))
                                (*checkStrength)[k] = false;
                            graphPtr->CalculateDistanceFromTo(subgraphPtr->vertices[i]->name, checkStrength);
                            if(graphPtr->vertices[subgraphPtr->vertices[i]->name]->distanceTo[subgraphPtr->vertices[j]->name] <= kParam){
                                if(numStrengthenedLazyCut == 0) aveNumTermsNegStrengthened = minimalSeparatorPtr->size();
                                else aveNumTermsNegStrengthened = (aveNumTermsNegStrengthened*numStrengthenedLazyCut + minimalSeparatorPtr->size())/(numStrengthenedLazyCut + 1);
                                numStrengthenedLazyCut++;
                            }

                            delete checkStrength;
                            delete minimalSeparatorPtr;
                            delete curSolCPrime;
                            delete toBeDeleted;

                        }else if(methodParam == "PGRAPHDELETEITERATE"){

                            vector<int>* curSolCPrime = new vector<int>();
                            vector<int>* toBeDeleted = new vector<int>();
                            for(int k = 0; k < curSolC->size(); k++){
                                int kVertex = (*curSolC)[k];
                                bool skipFlag = false;
                                for(auto gPtr : graphCollection){

                                    bool temp1 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[i]->name);
                                    bool temp2 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[j]->name);

                                    if(temp1 && temp2) continue;
                                    else{
                                        skipFlag = true;
                                        break;
                                    }
                                }

                                if(skipFlag == false) curSolCPrime->push_back(kVertex);
                                else toBeDeleted->push_back(kVertex);
                            }

                            vector<int>* minimalSeparatorPtr;
                            graph* sgPtr;
                            graph* gPtr;
                            int minimalSeparatorGraphIndex = -1;
                            for(int k = 0; k < graphCollection.size(); k++){

                                gPtr = graphCollection[k];

                                if(gPtr == graphPtr) sgPtr = subgraphPtr;
                                else sgPtr = gPtr->CreateSubgraph(curSol);

                                sgPtr->CalculateDistanceFrom(i);

                                if(sgPtr->vertices[i]->distanceTo[j] > kParam){

                                    vector<int>* tempSeparatorPtr = MINIMALIZE(gPtr, sgPtr->vertices[i]->name, sgPtr->vertices[j]->name, kParam, toBeDeleted, curSolCPrime);

                                    if(minimalSeparatorGraphIndex == -1){
                                        minimalSeparatorPtr = tempSeparatorPtr;
                                        minimalSeparatorGraphIndex = k;
                                    }else if(tempSeparatorPtr->size() < minimalSeparatorPtr->size()){
                                        delete minimalSeparatorPtr;
                                        minimalSeparatorPtr = tempSeparatorPtr;
                                        minimalSeparatorGraphIndex = k;
                                    }else delete tempSeparatorPtr;


                                }

                                if(sgPtr != subgraphPtr) delete sgPtr; // when sgPtr == subgraphPtr, we don't delete sgPtr here to avoid double delete

                            }

                            GRBLinExpr EXPR = X[subgraphPtr->vertices[i]->name] + X[subgraphPtr->vertices[j]->name];
                            for(auto k : (*minimalSeparatorPtr))
                                EXPR -= X[k];
                            addLazy(EXPR <= 1);

                            if(numLazyCut == 0) aveNumTermsNeg = minimalSeparatorPtr->size();
                            else aveNumTermsNeg = (aveNumTermsNeg*numLazyCut + minimalSeparatorPtr->size())/(numLazyCut + 1);

                            vector<bool>* checkStrength = new vector<bool>(graphPtr->n, true);
                            for(auto k : (*minimalSeparatorPtr))
                                (*checkStrength)[k] = false;

                            gPtr = graphCollection[minimalSeparatorGraphIndex];
                            gPtr->CalculateDistanceFromTo(subgraphPtr->vertices[i]->name, checkStrength);
                            if(gPtr->vertices[subgraphPtr->vertices[i]->name]->distanceTo[subgraphPtr->vertices[j]->name] <= kParam){
                                if(numStrengthenedLazyCut == 0) aveNumTermsNegStrengthened = minimalSeparatorPtr->size();
                                else aveNumTermsNegStrengthened = (aveNumTermsNegStrengthened*numStrengthenedLazyCut + minimalSeparatorPtr->size())/(numStrengthenedLazyCut + 1);
                                numStrengthenedLazyCut++;
                            }

                            delete checkStrength;
                            delete minimalSeparatorPtr;
                            delete curSolCPrime;
                            delete toBeDeleted;

                        }

                        numLazyCut++;
                        delete subgraphPtr;
                        goto theEnd;
                    }
                }
            }
            delete subgraphPtr;
        }

        theEnd:
        delete curSol;
        delete curSolC;

        callbackTime += GetWallTime() - timeBegin;

    }catch(GRBException e){
        cout << "ERROR NUMBER DURING : " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }catch(...){
        cout << "ERROR DURING CALLBACK" << endl;
    }
}



