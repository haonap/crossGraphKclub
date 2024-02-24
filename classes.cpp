//
// Created by Hao Pan on 9/20/21.
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

vector<int> graph::FindCommonNeighbors(int u, int v){
    vector<int> comNeighbors;
    int i = 0, j = 0;
    while(i < vertices[u]->neighbors.size() && j < vertices[v]->neighbors.size()){
        if(vertices[u]->neighbors[i] == vertices[v]->neighbors[j]){
            comNeighbors.push_back(vertices[u]->neighbors[i]);
            i++;
            j++;
        }else if (vertices[u]->neighbors[i] > vertices[v]->neighbors[j]){
            j++;
        }else{
            i++;
        }
    }
    return comNeighbors;
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
        int nbr_size = vertices[u]->neighbors.size();
        for(int i = 0; i < nbr_size; i++){
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

// remove edge uv from graph: first make sure if uv is an edge
void graph::DeleteEdge(int u, int v){
    vector<int>::iterator it = find(vertices[u]->neighbors.begin(), vertices[u]->neighbors.end(), v); //find location of vertex v
    //first check if uv is an edge; if not, return
    if (it==vertices[u]->neighbors.end())
        return;// v is NOT adjacent to u

    //remove v from u's neighbor list
    vertices[u]->neighbors.erase(it);
    vertices[u]->degree--;

    //remove u from v's neighbor list
   it = find(vertices[v]->neighbors.begin(), vertices[v]->neighbors.end(), u); //find location of vertex u
   //double check adjacency
    //double check if u is adjacent to v; if not, return
    if (it==vertices[v]->neighbors.end()) {
        cout<<"Caution!!! double check the graph!"<<endl;
        return;// v is NOT adjacent to u
    }
    vertices[v]->neighbors.erase(it);
    vertices[v]->degree--;

    m--; //edge decrease 1

    //update max degree and max maxDegreeVertexIndex
    maxDegree = -1;
    for (int i = 0; i < n; i++) {
        if (vertices[i]->degree > maxDegree){
            maxDegree = vertices[i]->degree;
            maxDegreeVertexIndex = i;
        }
    }
}

//graph remove nodes C: if node u in removedC, then removedC[u] = true;
// the function will reset the node degree =0 and removed it from its neighbors
void graph::GraphRemoveNodes(const vector<bool> &removedC){
    //reset removed node's degree
    for (int i =0; i < n; i++) {
        if (removedC[i]){
            //if vertex i in removedC, it means it will be removed: we reset its degree = 0
            vertices[i]->degree = 0;
            //remove i from all its neighbors
            for (int j = 0; j <vertices[i]->neighbors.size(); j++) {
                int i_Neighbor = vertices[i]->neighbors[j];// i's neighbor
                if (vertices[i_Neighbor]->degree > 0){
                    vertices[i_Neighbor]->degree--;
                }
            }
        }
    }
    //update each vertex's neighbor list
    vector<int> neighbors;
    int maxDeg = 0;
    int maxDegNode = -1;
    int edgeCount = 0;
    for(int i = 0; i < n; i++){
        neighbors.clear();
        if(vertices[i]->degree > 0){
            for(int j = 0; j < vertices[i]->neighbors.size(); j++){
                int pNeighbor = vertices[i]->neighbors[j];
                if(vertices[pNeighbor]->degree > 0){
                    neighbors.push_back(pNeighbor);
                }
            }
            vertices[i]->neighbors = neighbors;
            edgeCount += (int)neighbors.size();
            if(vertices[i]->degree > maxDeg){
                maxDeg = vertices[i]->degree;
                maxDegNode = i;
            }

        }else{
            vertices[i]->neighbors.clear();
        }
    }
    m = edgeCount/2;
    if(maxDeg > 0){
        maxDegree = maxDeg;
        maxDegreeVertexIndex = maxDegNode;
    }else{
        maxDegree = 0;
        maxDegreeVertexIndex = 0;
    }
}

// initialize graph, set #vertex = n, and all others 0 or empty
void graph::initialize(){
    for (int i =0; i < n; i++) {
        vertices[i]->neighbors.clear();
        vertices[i]->degree = 0;
        vertices[i]->kNeighbors.clear();
        vertices[i]->adjacency.clear();
        vertices[i]->commonNeighbors.clear();
        vertices[i]->numCommonNeighbors.size();
    }
    m = 0;
    maxDegree = 0;
    maxDegreeVertexIndex = 0;
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

    int maxDeg = 0;
    int maxDegreeVertexIndex = -1;
    for(int i = 0; i < graphPtr->n; i++){
        if(graphPtr->vertices[i]->degree > maxDeg){
            maxDeg = (int)graphPtr->vertices[i]->degree;
            maxDegreeVertexIndex = i;
        }
    }
    graphPtr->maxDegree = maxDeg;
    graphPtr->maxDegreeVertexIndex = maxDegreeVertexIndex;

    return graphPtr;
}

// find all vertices connected to a vertex
vector<int> graph::BFS(int s){
    vector<int> tempVec;
    vector<int> distance(n, n);
    queue<int> Q;
    distance[s] = 0;
    Q.push(s);
    int u, v;
    int goFlag = 1;

    while(Q.size() && goFlag){
        u = Q.front();
        Q.pop();
        for(int i = 0; i < vertices[u]->neighbors.size(); i++){
            v = vertices[u]->neighbors[i];
            if(distance[v] > distance[u] + 1){
                distance[v] = distance[u] + 1;
                Q.push(v);
            }
        }
    }

    for(int i = 0; i < distance.size(); i++){
        if(distance[i] < n && i != s){
            tempVec.push_back(i);
        }
    }
    return tempVec;
}


// calculate distance from source to all other nodes in a graph
void graph::CalculateDistanceFrom(int source){
    vector<bool>* V = new vector<bool>(n, true);
    CalculateDistanceFromTo(source, V);
}

//calculate distance from source to other nodes in subgraph G[S]
void graph::CalculateDistanceFromTo(int source, vector<bool>* S) {
    queue<int>* Q = new queue<int>();
    vertices[source]->distanceTo.clear();
    vertices[source]->distanceTo.resize(n, n*1000); //in case graph is not connected, we should assign a big number.
    vertices[source]->distanceTo[source] = 0;
    Q->push(source);
    int u, v;

    while(Q->size()){
        u = Q->front();
        Q->pop();
        for(int i = 0; i < vertices[u]->degree; i++){
            v = vertices[u]->neighbors[i];
            if((*S)[v] == false) continue; //if v is NOT in G[S], skip this v
            if(vertices[source]->distanceTo[v] > vertices[source]->distanceTo[u] + 1){
                vertices[source]->distanceTo[v] = vertices[source]->distanceTo[u] + 1;
                Q->push(v);
            }
        }
    }

    delete Q;
}


pGraph_callback::pGraph_callback(string methodParam, GRBVar* X, vector<graph*>& graphCollection, int kParam, int& numStrengthenedLazyCut, int& numLazyCut, int& numCallback, double& callbackTime, double& countSLC_Time,double& aveNumTermsNeg, double& aveNumTermsNegStrengthened, vector<GRBModel>& SEP_MODELS, vector<GRBVar*>& Z, vector<GRBVar*>& W, double epsilon, int& independentSetInequalityCount, int& independentSetInequalityCountViolateEpsilon) : methodParam(methodParam), X(X), graphCollection(graphCollection), kParam(kParam), numStrengthenedLazyCut(numStrengthenedLazyCut), numLazyCut(numLazyCut), numCallback(numCallback), callbackTime(callbackTime),countSLC_Time(countSLC_Time), aveNumTermsNeg(aveNumTermsNeg), aveNumTermsNegStrengthened(aveNumTermsNegStrengthened), SEP_MODELS(SEP_MODELS), Z(Z), W(W), epsilon(epsilon), independentSetInequalityCount(independentSetInequalityCount), independentSetInequalityCountViolateEpsilon(independentSetInequalityCountViolateEpsilon){}


void pGraph_callback::callback(){
    try {

        if (where != GRB_CB_MIPSOL) return;
        numCallback++;
        double timeBegin = GetWallTime();
        int nParam = graphCollection[0]->n;
        int pParam = graphCollection.size(); //number of colllections
        double *x = getSolution(X, nParam);
        vector<int> *curSol = new vector<int>();
        vector<int> *curSolC = new vector<int>();
        for (int i = 0; i < nParam; i++) {
            if (x[i] > 0.5) curSol->push_back(i);
            else curSolC->push_back(i);
        }

        if (curSol->empty()){
            delete curSol;
            delete curSolC;
            return;
        }

        // use PPCF + Independent set inequality: regardless of IS adding cut, we still need PPCF cut
        if (methodParam == "PPCF_IS") {
            //add Independent set cut starting here
            for (int i = 0; i < SEP_MODELS.size(); i++) {
                GRBModel cur_model = SEP_MODELS[i];

                GRBVar *temp = cur_model.getVars();

                for (int j = 0; j < nParam * 2; j++) {
                    string group = temp[j].get(GRB_StringAttr_VarName).substr(0, 1);
                    string index = temp[j].get(GRB_StringAttr_VarName).substr(1);

                    if (find(curSol->begin(), curSol->end(), stoi(index)) != curSol->end()) {
                        if (group == "W")
                            temp[j].set(GRB_DoubleAttr_Obj, -1.0);
                        if (group == "Z")
                            temp[j].set(GRB_DoubleAttr_Obj, 1.0);
                    }
                }

                cur_model.optimize();
                if (cur_model.get(GRB_DoubleAttr_ObjVal) > 1) independentSetInequalityCount++;
                if (cur_model.get(GRB_DoubleAttr_ObjVal) >= 1 + epsilon) {
                    //cout<<"IS cut is added"<<endl;
                    GRBLinExpr EXPR = 0;
                    if (cur_model.get(GRB_IntAttr_SolCount)) {
                        for (int j = 0; j < nParam * 2; j++) {
                            string group = temp[j].get(GRB_StringAttr_VarName).substr(0, 1);
                            string index = temp[j].get(GRB_StringAttr_VarName).substr(1);
                            if (temp[j].get(GRB_DoubleAttr_X) > 0) {
                                if (group == "Z") {
                                    EXPR += X[stoi(index)] * temp[j].get(GRB_DoubleAttr_X);
                                }
                                if (group == "W") {
                                    EXPR -= X[stoi(index)] * temp[j].get(GRB_DoubleAttr_X);
                                }
                            }
                        }
                    }
                    addLazy(EXPR <= 1);
                    independentSetInequalityCountViolateEpsilon++;
                    break;
                }
            }
            //add Independent set cut ending here
        }//end if

        //code updated by Yl
        int n_curSol = curSol->size();
        vector<graph*> subGraphCollection; //store sub graphs induced by curSol
        //calculate distance from i to others in the subgraphs
        for (auto graphPtr1: graphCollection) {
            graph * subgraphPtr1 = graphPtr1->CreateSubgraph(curSol);
            for (int i = 0; i < n_curSol; i++) {
                subgraphPtr1->CalculateDistanceFrom(i);
            }
            subGraphCollection.push_back(subgraphPtr1);
        }


        graph *subgraphPtr;
        graph *graphPtr;
        //add lazy cut
        for (int i = 0; i < n_curSol; i++) {
            int u = (*curSol)[i];
            for (int j = i + 1; j < n_curSol; j++) {
                int v = (*curSol)[j];
                bool isPeeled = false;//for each pair of uv, only needs to peel graphs once
                vector<int> *curSolCPrime = new vector<int>();
                vector<int> *toBeDeleted = new vector<int>();
                for (int k1 = 0; k1 < pParam; k1++) {
                    graphPtr = graphCollection[k1];
                    subgraphPtr = subGraphCollection[k1];
                    //check if there is need to add cut
                    if (subgraphPtr->vertices[i]->distanceTo[j] > kParam){
                        //CCF cut
                        if (methodParam == "BASE") {
                            vector<int> *minimalSeparatorPtr = MINIMALIZE(graphPtr, subgraphPtr->vertices[i]->name,
                                                                          subgraphPtr->vertices[j]->name, kParam,curSolC);
                            GRBLinExpr EXPR = X[subgraphPtr->vertices[i]->name] + X[subgraphPtr->vertices[j]->name];
                            for (int k = 0; k < minimalSeparatorPtr->size(); k++)
                                EXPR -= X[(*minimalSeparatorPtr)[k]];
                            addLazy(EXPR <= 1);

                            if (numLazyCut == 0) aveNumTermsNeg = minimalSeparatorPtr->size();
                            else
                                aveNumTermsNeg =(aveNumTermsNeg * numLazyCut + minimalSeparatorPtr->size()) / (numLazyCut + 1);
                            delete minimalSeparatorPtr;
                        //end if = 'Base'
                        } else{
                            //This else condition is for the methodParam == "PPCF" or " PPCF_IS"
                            //when method PPCF or PPCF_IS, we need to add PPCF cut

                            //only need to peel once for each pair of u and v
                            if (!isPeeled){
                                isPeeled = true;
                                //only need to peel once for each pair of u and v
                                for (int k = 0; k < curSolC->size(); k++) {
                                    int kVertex = (*curSolC)[k];
                                    bool skipFlag = false;
                                    for (auto gPtr: graphCollection) {
                                        bool temp1 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[i]->name);
                                        bool temp2 = gPtr->IsKAdjacent(kVertex, subgraphPtr->vertices[j]->name);
                                        if (temp1 && temp2) continue;
                                        else {
                                            skipFlag = true;
                                            break;
                                        }
                                    }
                                    if (skipFlag == false) curSolCPrime->push_back(kVertex);
                                    else toBeDeleted->push_back(kVertex);
                                }

                                // recursively delete vertices which are too far away from i or j (actually only removing edges in implementation)
                                vector<bool> *toBeKeptBool = new vector<bool>(graphPtr->n, true);
                                for (auto q: *toBeDeleted)
                                    (*toBeKeptBool)[q] = false;

                                int numVertexDeleted;
                                do {
                                    numVertexDeleted = 0;
                                    for (auto gPtr: graphCollection) {
                                        gPtr->CalculateDistanceFromTo(subgraphPtr->vertices[i]->name, toBeKeptBool);
                                        gPtr->CalculateDistanceFromTo(subgraphPtr->vertices[j]->name, toBeKeptBool);
                                        for (int p = 0; p < gPtr->n; p++) {
                                            if (p == subgraphPtr->vertices[i]->name ||
                                                p == subgraphPtr->vertices[j]->name)
                                                continue;
                                            if (gPtr->vertices[subgraphPtr->vertices[i]->name]->distanceTo[p] > kParam) {
                                                if ((*toBeKeptBool)[p] == true) {
                                                    (*toBeKeptBool)[p] = false;
                                                    numVertexDeleted++;
                                                }
                                            }

                                            if (gPtr->vertices[subgraphPtr->vertices[j]->name]->distanceTo[p] > kParam) {
                                                if ((*toBeKeptBool)[p] == true) {
                                                    (*toBeKeptBool)[p] = false;
                                                    numVertexDeleted++;
                                                }
                                            }

                                        }
                                    }
                                } while (numVertexDeleted > 0);
                                //end iterative peeling

                                curSolCPrime->clear();
                                toBeDeleted->clear();
                                for (int k = 0; k < toBeKeptBool->size(); k++) {
                                    if ((*toBeKeptBool)[k] == false) toBeDeleted->push_back(k);
                                    else if (find(curSol->begin(), curSol->end(), k) == curSol->end())
                                        curSolCPrime->push_back(k);
                                }
                                delete toBeKeptBool;
                            }//end if for checking if is_Peel

                            // find PPCF and add lazy cut
                            vector<int> *minimalSeparatorPtr = MINIMALIZE(graphPtr, subgraphPtr->vertices[i]->name,
                                                                          subgraphPtr->vertices[j]->name, kParam, toBeDeleted, curSolCPrime);

                            GRBLinExpr EXPR = X[subgraphPtr->vertices[i]->name] + X[subgraphPtr->vertices[j]->name];
                            for (auto k: (*minimalSeparatorPtr))
                                EXPR -= X[k];
                            addLazy(EXPR <= 1);
                            // calculate aveNumTermsNeg
                            if (numLazyCut == 0)
                                aveNumTermsNeg = minimalSeparatorPtr->size();
                            else
                                aveNumTermsNeg = (aveNumTermsNeg * numLazyCut + minimalSeparatorPtr->size()) / (numLazyCut + 1);

                            //check if PPCF is strengthened cut, i.e., check if scuh PPCF cut is NOT CCF aross ALL original graph collection
                            double SLCtimeBegin = GetWallTime(); //  time for counting number of SLC
                            vector<bool> *checkStrength = new vector<bool>(graphPtr->n, true);
                            for (auto k: (*minimalSeparatorPtr))
                                (*checkStrength)[k] = false;

                            graph *gPtr;
                            bool isSLC = true; //used to mark if such PPCF is strengthened cut
                            for (int k = 0; k < graphCollection.size(); k++) {
                                gPtr = graphCollection[k];
                                gPtr->CalculateDistanceFromTo(subgraphPtr->vertices[i]->name, checkStrength);
                                if (gPtr->vertices[subgraphPtr->vertices[i]->name]->distanceTo[subgraphPtr->vertices[j]->name] >
                                    kParam) {
                                    isSLC = false;
                                    break;
                                }
                            }
                            countSLC_Time += GetWallTime() - SLCtimeBegin; //time for counting number of SLC
                            // if isSLC is true, it means that such PPCF is strengthened
                            if (isSLC){
                                numStrengthenedLazyCut++;
                                if (numStrengthenedLazyCut == 0)
                                    aveNumTermsNegStrengthened = minimalSeparatorPtr->size();
                                else
                                    aveNumTermsNegStrengthened = (aveNumTermsNegStrengthened * numStrengthenedLazyCut +
                                                                  minimalSeparatorPtr->size()) /(numStrengthenedLazyCut + 1);
                            }

                            delete checkStrength;
                            delete minimalSeparatorPtr;

                        }//end else
                        numLazyCut++;
                        //uncomment below so that only one cut per callback; otherwise comment below to run the version: add more cuts per callback
                        //delete curSolCPrime;
                        //delete toBeDeleted;
                        //goto theEnd;
                    } //end if >kParam
                }//end for loop with k1< pParam
                delete curSolCPrime;
                delete toBeDeleted;
            }//end for loop with j <n_curSol
        }//end for loop with i <n_curSol

        //theEnd:
        delete curSol;
        delete curSolC;

        //clear subgraph collection
        for(auto subgraphPtr1 : subGraphCollection)
            delete subgraphPtr1;
        subGraphCollection.clear();

        //cout<<"finished one callback and deleted subGraphCollection"<<endl;

        callbackTime += GetWallTime() - timeBegin;

    }catch(GRBException e){
        cout << "ERROR NUMBER DURING : " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }catch(...){
        cout << "ERROR DURING CALLBACK" << endl;
    }
}

vector<int> graph::FindDegeneracyOrdering(vector<int>& rightDegree){
    int degeneracy = 0;
    rightDegree.resize(n);

    for(int i = 0; i < n; i++)
        rightDegree[i] = vertices[i]->degree;

    vector<int> bin(maxDegree + 1, 0);
    for(int i = 0; i < n; i++) bin[rightDegree[i]]++;
    int start = 0;
    for(int d = 0; d <= maxDegree; d++){
        int num = bin[d];
        bin[d] = start;
        start += num;
    }

    vector<int> pos(n);
    vector<int> vert(n);
    for(int i = 0; i < n; i++){
        pos[i] = bin[rightDegree[i]];
        vert[pos[i]] = i;
        bin[rightDegree[i]]++;
    }

    for(int d = maxDegree; d >=1; d--) bin[d] = bin[d - 1];
    bin[0] = 0;

    for(int i = 0; i < n; i++){
        int minv = vert[i];
        bin[rightDegree[minv]]++;
        degeneracy = max(degeneracy, rightDegree[minv]);
        for(int j = 0; j < vertices[minv]->degree; j++){
            int u = vertices[minv]->neighbors[j];
            if(pos[u] > pos[minv]){
                if(rightDegree[u] == rightDegree[minv]){
                    int pw = bin[rightDegree[minv]];
                    int w = vert[pw];
                    if(u != w){
                        vert[pw] = u;
                        vert[pos[u]] = w;
                        pos[w] = pos[u];
                        pos[u] = pw;
                    }
                    bin[rightDegree[minv] - 1] = pos[minv] + 1;
                    bin[rightDegree[u]]++;
                    rightDegree[u]--;
                }else{
                    int pw = bin[rightDegree[u]];
                    int w =vert[pw];

                    if(u != w){
                        vert[pw] = u;
                        vert[pos[u]] = w;
                        pos[w] = pos[u];
                        pos[u] = pw;
                    }
                    bin[rightDegree[u]]++;
                    rightDegree[u]--;
                }
            }
        }
    }
    return vert;
}

vector<int> graph::FindHeuristicClique(vector<int>& degeneracyOrder, vector<int>& rightDegree) {
    vector<int> clique;
    for(int i = 0; i < n && clique.empty(); i++){
        int v = degeneracyOrder[i];
        if(rightDegree[v] == n - i - 1){
            clique.resize(n - i);
            for(int j = i; j < n; j++) clique[j - i] = degeneracyOrder[j];
            sort(clique.begin(), clique.end());
        }
    }
    return clique;
}

vector<int> graph::ShortestPathsUnweighted(int origin, vector<bool>& S) {
    vector<int> D;
    D.push_back(origin);
    return MultiSourceShortestPaths(D, S);
}

vector<int> graph::MultiSourceShortestPaths(vector<int>& D, vector<bool>& S) {
    long u,v;
    vector<int> dist(n, n);
    vector<bool> reached(n, false);
    vector<int> children, parents;
    bool status = false;

    for(int i = 0; i < D.size(); i++){
        if(S[D[i]]) status = true;
        else continue;

        children.push_back(D[i]);
        dist[D[i]] = 0;
        reached[D[i]] = true;
    }

    if(!status) return dist;

    for(int d = 1; !children.empty(); d++){
        parents = children;
        children.clear();
        for(int i = 0; i < parents.size(); i++){
            u = parents[i];
            for(int j = 0; j < vertices[u]->degree; j++){
                v = vertices[u]->neighbors[j];
                if(!reached[v] && S[v]){
                    reached[v] = true;
                    dist[v] = d;
                    children.push_back(v);
                }
            }
        }
    }
    return dist;
}



