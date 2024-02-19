//
// Created by Hao Pan on 9/20/21.
//

#include "functions.h"
#include "classes.h"

//Color::Modifier green(Color::FG_GREEN);
//Color::Modifier def(Color::FG_DEFAULT);


vector<string>* ReadTasks(){

    vector<string>* tasks = new vector<string>();
    string temp;
    ifstream fin("tasksCrossGraph.txt");
    while(getline(fin, temp))
        tasks->push_back(temp);
    return tasks;
}

vector<graph*> graphSequence;
vector<graph*> graphCollection;
int pParam, kParam, tauParam, TParam, firstGraph, lastGraph;
long num_edges_collection_before_peel = 0; //before peeling, the total number of edges across all grpahs in the collection
long num_edges_collection_after_peel = 0; //after peeling, the total number of edges across all grpahs in the collection
int num_vertex_before_peel = 0; //before peeling, the total number of vertices
int num_vertex_after_peel = 0; //after peeling, the total number of vertices
string instanceName, methodParam;
double epsilonParam;
void ReadInput(string task){

    stringstream ss(task);
    int u, v;
    string temp;

    ss >> pParam >> firstGraph >> lastGraph >> instanceName >> kParam >> methodParam >> epsilonParam;
    //cout << green << "\nREAD IN INPUT..." << def << endl;
    cout << "\nREAD IN INPUT..." << endl;
    cout <<  "p = " << pParam << ", k = " << kParam << endl;

    //check if the condition for using PPCF_IS
    if (methodParam =="PPCF_IS"){
        if (kParam != 2){
            cout<<"PPCF_IS method is ONLY works for k=2!!! Please use BASE or PPCF when k is not equal to 2!"<<endl;
            exit(0);
        }
    }

    cout << "INSTANCE : " << instanceName << " " << firstGraph << "-" << lastGraph << endl;
    cout << "EPSILON FOR IS SEPARATION : " << epsilonParam << endl;
    cout<<"The method for solving instance-"<<instanceName << " " << firstGraph << "-" << lastGraph <<" is: "<<methodParam<<endl;

    ifstream fin;

    //initialization
    num_edges_collection_before_peel = 0; //before peeling, the total number of edges across all grpahs in the collection
    num_edges_collection_after_peel = 0; //after peeling, the total number of edges across all grpahs in the collection
    num_vertex_before_peel = 0; //before peeling, the total number of vertices
    num_vertex_after_peel = 0; //after peeling, the total number of vertices

    for(int i = 0; i < pParam; i++){
        graph* graphPtr = new graph();
        string graphName = "./instances/" + instanceName + "/" + instanceName + "_" + itos_c(i + firstGraph) + ".txt";
        fin.open(graphName, ifstream::in);
        fin >> temp >> graphPtr->n >> temp >> graphPtr->m;

        num_vertex_before_peel = graphPtr->n;// store total number of edges before peeling: every graph has the same set of nodes
        num_edges_collection_before_peel += graphPtr->m;// store total number of edges before peeling

        graphPtr->ind = i;
        for(int j = 0; j < graphPtr->n; j++){
            vertex* vertexPtr = new vertex();
            vertexPtr->ind = j;
            graphPtr->vertices.push_back(vertexPtr);
        }

        while(fin >> temp){
            if(temp == "e"){
                fin >> u >> v;
                graphPtr->vertices[u - 1]->neighbors.push_back(v - 1);
                graphPtr->vertices[v - 1]->neighbors.push_back(u - 1);
            }
        }

        graphPtr->maxDegree = 0;
        for(int j = 0; j < graphPtr->n; j++){
            graphPtr->vertices[j]->degree = graphPtr->vertices[j]->neighbors.size();
            if(graphPtr->vertices[j]->degree > graphPtr->maxDegree){
                graphPtr->maxDegree = graphPtr->vertices[j]->degree;
                graphPtr->maxDegreeVertexIndex = j;
            }
        }
        fin.close();
        graphCollection.push_back(graphPtr);
    }
}

void ReadInput(int first, int last){

    int u, v;
    string temp;

    firstGraph = first;
    lastGraph = last;

    cout << "\nREAD IN INPUT..." << endl;
    cout <<  "p = " << pParam << ", k = " << kParam << endl;
    cout << "INSTANCE : " << instanceName << " " << firstGraph << "-" << lastGraph << endl;

    ifstream fin;

    for(int i = 0; i < pParam; i++){
        graph* graphPtr = new graph();
        string graphName = "./instances/" + instanceName + "/" + instanceName + "_" + itos_c(i + firstGraph) + ".txt";
        fin.open(graphName, ifstream::in);
        fin >> temp >> graphPtr->n >> temp >> graphPtr->m;

        graphPtr->ind = i;
        for(int j = 0; j < graphPtr->n; j++){
            vertex* vertexPtr = new vertex();
            vertexPtr->ind = j;
            graphPtr->vertices.push_back(vertexPtr);
        }

        while(fin >> temp){
            if(temp == "e"){
                fin >> u >> v;
                graphPtr->vertices[u - 1]->neighbors.push_back(v - 1);
                graphPtr->vertices[v - 1]->neighbors.push_back(u - 1);
            }
        }

        graphPtr->maxDegree = 0;
        for(int j = 0; j < graphPtr->n; j++){
            graphPtr->vertices[j]->degree = graphPtr->vertices[j]->neighbors.size();
            if(graphPtr->vertices[j]->degree > graphPtr->maxDegree){
                graphPtr->maxDegree = graphPtr->vertices[j]->degree;
                graphPtr->maxDegreeVertexIndex = j;
            }
        }
        fin.close();
        graphCollection.push_back(graphPtr);
    }
}

void ReadInputKClubSig(string task){

    stringstream ss(task);

    ss >> instanceName >> tauParam >> kParam >> TParam >> methodParam;
    //cout << green << "\nREAD IN INPUT..." << def << endl;
    cout << "\nREAD IN INPUT..." << endl;
    cout <<  "tau = " << tauParam << ", k = " << kParam << ", T = " << TParam << endl;
    cout << "INSTANCE : " << instanceName << endl;
}


// get a heuristic p-graph k-club, actually a k-club on the intersection graph of all graphs in the collection
vector<int> GetHeuristicSol(){
    vector<int> HeuSol = {};

    graph* powerGraphPtr = GetIntersectionGraphOfPowerGraphs();
    vector<int> rd; // right-degree of degeneracy ordering
    vector<int> degeneracyOrdering = powerGraphPtr->FindDegeneracyOrdering(rd);
    vector<int> kclique = powerGraphPtr->FindHeuristicClique(degeneracyOrdering, rd);

    vector<bool> heuristicSolution = boolify(kclique, powerGraphPtr->n);
    graph* graphPtr = GetIntersectionGraph();
    HeuSol = Drop(graphPtr, heuristicSolution);

//    graph* graphPtr = GetIntersectionGraph();
//
//    HeuSol = graphPtr->vertices[graphPtr->maxDegreeVertexIndex]->neighbors;
//    HeuSol.push_back(graphPtr->maxDegreeVertexIndex);
//    sort(HeuSol.begin(), HeuSol.end());

    delete powerGraphPtr;
    delete graphPtr;
    return HeuSol;
}

vector<int> Drop(graph* graphPtr, vector<bool>& W){
    vector<int> near(graphPtr->n, 0);

    int Wsize = count(W.begin(), W.end(), true);

    for(int size = Wsize; size >= 0; size--){
        for(int i = 0; i < graphPtr->n; i++){
            if(!W[i]) continue;
            near[i] = KNeighborhoodSize(graphPtr, W, i);
        }

        int w;
        int smallestNearby = size;
        for(int i = 0; i < graphPtr->n; i++){
            if(!W[i]) continue;
            if(near[i] < smallestNearby){
                w = i;
                smallestNearby = near[i];
            }
        }

        if(smallestNearby == size) break;
        W[w] = false;
    }

    vector<int> Wvector;
    for(int i = 0; i < graphPtr->n; i++){
        if(W[i]){
            Wvector.push_back(i);
        }
    }

    return Wvector;
}

int KNeighborhoodSize(graph* graphPtr, vector<bool>& W, int v){
    if(!W[v]){
        cerr << "\n ERROR: K-neighborhoodSize. You are calculating distances across W nodes starting from some vertex v which does not belong to W.";
        return 0;
    }
    vector<int> dist = graphPtr->ShortestPathsUnweighted(v, W);
    int nsize = 0;
    for(int i = 0; i < graphPtr->n; i++) if(dist[i] <= kParam) nsize++;
    return nsize;
}



graph* GetIntersectionGraph(){

    if(pParam == 1){
        graph* graphPtr = new graph();
        *graphPtr = *(graphCollection[0]);
        return graphPtr;
    }else{
        graph* graphPtr = new graph();
        graphPtr->n = graphCollection[0]->n;
        for(int i = 0; i < graphPtr->n; i++){
            vertex* vertexPtr = new vertex();
            vertexPtr->ind = i;
            graphPtr->vertices.push_back(vertexPtr);
        }

        int edgeCount = 0;
        for(int i = 0; i < graphPtr->n; i++){
            for(int j = 0; j < graphCollection[0]->vertices[i]->neighbors.size(); j++){
                int pNeighbor = graphCollection[0]->vertices[i]->neighbors[j];
                if(pNeighbor > i){
                    int graphCount = 1;
                    for(int t =  1; t < graphCollection.size(); t++){
                        if(graphCollection[t]->IsAdjacent(i, pNeighbor)){
                            graphCount++;
                        }else{
                            graphCount = 0;
                            break;
                        }
                    }
                    if(graphCount > 0){
                        graphPtr->vertices[i]->neighbors.push_back(pNeighbor);
                        graphPtr->vertices[pNeighbor]->neighbors.push_back(i);
                        edgeCount++;
                    }
                }
            }
        }

        graphPtr->m = edgeCount;

        int maxDeg = 0;
        int maxDegreeVertexIndex = -1;
        for(int i = 0; i < graphPtr->n; i++){
            graphPtr->vertices[i]->degree = (int)graphPtr->vertices[i]->neighbors.size();
            if(graphPtr->vertices[i]->degree > maxDeg){
                maxDeg = (int)graphPtr->vertices[i]->degree;
                maxDegreeVertexIndex = i;
            }
        }
        graphPtr->maxDegree = maxDeg;
        graphPtr->maxDegreeVertexIndex = maxDegreeVertexIndex;
        return graphPtr;
    }
}

vector<int> GetPersistentKClubSig(int windowHead, int* peelPtr, vector<int>* flagOptimalPtr){
    double timeBegin = GetWallTime();
//    int upperBound = (int)graphCollection[0]->n + 1;
    vector<int> solution;
    solution.clear();

    vector<int> HeuSol;
    HeuSol.clear();
    if(windowHead == 0) {
        HeuSol = GetHeuristicSol();
        if(*peelPtr < HeuSol.size()) *peelPtr = HeuSol.size();
    }


    graph* intersectionOfPowerGraphPtr;
    bool isChanged = true;
    int edge_after_peel = 0;
    double peelTimeBegin = GetWallTime();
    //Recursive peeling--Algo 2
    while(isChanged){
        isChanged = false;
        //find total number of edges before peel
        int edge_before_peel = 0;
        for (auto g: graphCollection) {
            edge_before_peel += g->m;
        }
        intersectionOfPowerGraphPtr = GetIntersectionGraphOfPowerGraphs(); //obtain the power intersection graph J ( Gk )
        CorePeel(intersectionOfPowerGraphPtr, *peelPtr);// core peeling
        CommunityPeel(intersectionOfPowerGraphPtr, *peelPtr);//Community Peel
        DisconnectGraphsByComponents(graphCollection, intersectionOfPowerGraphPtr);//Cross Edge Peel

        //find total number of edges after peel
        edge_after_peel = 0;
        for (auto g: graphCollection) {
            edge_after_peel += g->m;
        }

        if (edge_before_peel != edge_after_peel)
            isChanged = true;
    }
    double peelTimeEnd = GetWallTime();
    //As the optimal solutions only exist in components whose size >= *peelPtr
    vector<vector<int>> componentsAfterPeeling = GetComponents(intersectionOfPowerGraphPtr, *peelPtr);
    // get vertex set after peeling
    vector<int> verticesAfterPeeling;
    for(int i = 0; i < componentsAfterPeeling.size(); i++){
        for(int j = 0; j < componentsAfterPeeling[i].size(); j++){
            verticesAfterPeeling.push_back(componentsAfterPeeling[i][j]);
        }
    }
    //sort
    sort(verticesAfterPeeling.begin(),verticesAfterPeeling.end());
    // get intersection of power graph after peeling
    graph* intersectionOfPowerGraphPtrAfterPeeling = intersectionOfPowerGraphPtr->CreateSubgraph(&verticesAfterPeeling);
    vector<graph*> graphCollectionAfterPeeling;

    //get number of vertices and edges after peeling
    num_vertex_after_peel = verticesAfterPeeling.size();
    num_edges_collection_after_peel = 0;
    graph *g;
    for(auto i : graphCollection){
        g = i->CreateSubgraph(&verticesAfterPeeling);
        num_edges_collection_after_peel += g->m;
        graphCollectionAfterPeeling.push_back(g);
    }

    //get components of power intersection graph
    componentsAfterPeeling = GetComponents(intersectionOfPowerGraphPtrAfterPeeling, 0);

    double timeMiddle = GetWallTime();

    if(componentsAfterPeeling.size() == 0) {
        cout << "NOTHING LEFT AFTER PREPROCESSING" << endl;
        return HeuSol;
    }

    try{

        for(auto graphPtr : graphCollectionAfterPeeling)
            graphPtr->FindKNeighbors(kParam);

        GRBEnv* ENV = 0;
        GRBVar* X = 0;
        GRBVar* Y = 0;

        ENV = new GRBEnv();
        GRBModel MODEL = GRBModel(*ENV);
        int nParam = graphCollectionAfterPeeling[0]->n;

        X = MODEL.addVars(nParam, GRB_BINARY);
        Y = MODEL.addVars(componentsAfterPeeling.size(), GRB_BINARY);
        MODEL.update();

        for(int i = 0; i < nParam; i++)
            X[i].set(GRB_DoubleAttr_Obj, 1.0);
        MODEL.update();

        GRBLinExpr sumY;
        for(int i = 0; i < componentsAfterPeeling.size(); i++){
            sumY += Y[i];
            for(int j = 0; j < componentsAfterPeeling[i].size(); j++){
                MODEL.addConstr(X[componentsAfterPeeling[i][j]] <= Y[i]);
            }
        }
        MODEL.addConstr(sumY == 1);

        for(int i = 0; i < componentsAfterPeeling.size(); i++){
            for(int u = 0; u < componentsAfterPeeling[i].size(); u++){
                int uVertex = componentsAfterPeeling[i][u];
                for(int v = u + 1; v < componentsAfterPeeling[i].size(); v++){
                    int vVertex = componentsAfterPeeling[i][v];
                    if(intersectionOfPowerGraphPtrAfterPeeling->IsAdjacent(uVertex,vVertex) == false)
                        MODEL.addConstr(X[uVertex] + X[vVertex] <= 1);
                }
            }
        }
        MODEL.update();
        /*
        MODEL.set(GRB_IntParam_MIPFocus, 0);
        MODEL.set(GRB_IntParam_Threads, 0);
        MODEL.set(GRB_IntParam_Method, -1);
        MODEL.set(GRB_IntParam_NodeMethod, 1);
         */
        MODEL.set(GRB_IntParam_Cuts, 0);

        MODEL.set(GRB_DoubleParam_TimeLimit, 3600);
        //Set termination gap limit; as needed; default is 1e-4
        //MODEL.set(GRB_DoubleParam_MIPGap, 1e-4);
        MODEL.set(GRB_StringParam_LogFile, "");
        MODEL.set(GRB_IntParam_OutputFlag, 1);
        MODEL.set(GRB_IntAttr_ModelSense, -1);
        MODEL.set(GRB_IntParam_LazyConstraints,1);
        int numStrengthenedLazyCut = 0, numLazyCut = 0, numCallback = 0, independentSetInequalityCount = 0, independentSetInequalityCountViolateEpsilon = 0;
        double callbackTime = 0, aveNumTermsNeg = 0, aveNumTermsNegStrengthened = 0;

        //if Method does not use PPCF_IS, they are empty set; this is to have a general callback function
        vector<GRBModel> SEP_MODELS;
        vector<GRBVar*> Z(pParam);
        vector<GRBVar*> W(pParam);

        pGraph_callback myCallback = pGraph_callback(methodParam, X, graphCollectionAfterPeeling, kParam, numStrengthenedLazyCut, numLazyCut, numCallback, callbackTime, aveNumTermsNeg, aveNumTermsNegStrengthened,SEP_MODELS, Z, W, 0.5, independentSetInequalityCount, independentSetInequalityCountViolateEpsilon);
        MODEL.setCallback(&myCallback);

        MODEL.optimize();

        double timeEnd = GetWallTime();

        //validate solution
        vector<int> *solutionPtr1 = new vector<int>();
        if(MODEL.get(GRB_IntAttr_SolCount)) {
            for (int i = 0; i < nParam; i++) {
                if (X[i].get(GRB_DoubleAttr_X) > 0.5){
                    solution.push_back(intersectionOfPowerGraphPtrAfterPeeling->vertices[i]->name);
                    solutionPtr1->push_back(intersectionOfPowerGraphPtrAfterPeeling->vertices[i]->name);
                }
            }
            //validate solution

            if(!ValidateSolution(solutionPtr1)){
                cout<<"The solution is invalid!!!! Caution!!!"<<endl;
                delete solutionPtr1;
            }
            delete solutionPtr1;

            cout << "OPTIMALITY STATUS: " << (MODEL.get(GRB_IntAttr_Status) == GRB_OPTIMAL) << endl;
            if(MODEL.get(GRB_IntAttr_Status) != GRB_OPTIMAL) (*flagOptimalPtr)[windowHead] = 0;
            if(*peelPtr < (int)solution.size()) *peelPtr = (int)solution.size();
            delete ENV;
            if((int)solution.size() > (int)HeuSol.size()) return solution;
            else return HeuSol;
        }

        delete ENV;

    }catch (GRBException e) {
        cout << "ERROR CODE = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "EXCEPTION DURING OPTIMIZATION" << endl;
    }

    (*flagOptimalPtr)[windowHead] = 0;
    return HeuSol;
}


void KClubSig(){
    cout << "\nSOLVE FOR A MAX " << tauParam << "-PERSISTENT " << kParam << "-CLUB SIG... " << endl;
    cout << "\nMETHOD: " << methodParam << endl;
    double st = GetWallTime();
    vector<int> bestKClub;
    bestKClub.clear();
    vector<int> flagOptimal(TParam - tauParam + 1, 1);
    int peel = 0;
    int bestWindow = 0;
    pParam = tauParam;

    for(int i = 0; i < TParam - tauParam + 1; i++) {

        cout << "\nIN WINDOW " << i + 1 << "..." << endl;

        ReadInput(i + 1, i + tauParam);

        vector<int> pKClub;
        pKClub.clear();

        pKClub = GetPersistentKClubSig(i, &peel, &flagOptimal);

        if (bestKClub.size() < pKClub.size()) {
            bestKClub = pKClub;
            bestWindow = i;
        }

        if (i > 0) {
            if (flagOptimal[i] == 0) {
                if (flagOptimal[i - 1] == 0) {
                    CleanUp();
                    break;
                }
            }
        }

        CleanUp();
    }

    cout << "\nBEST PERSISTENT KCLUB FOUND: ";
    int flagOpt = 1;//Check if instance is solved to optimal: flagOpt =1, optimal;
    for(int i = 0; i < flagOptimal.size(); i++) flagOpt *= flagOptimal[i];
    for(int i = 0; i < bestKClub.size(); i++) cout << bestKClub[i] + 1 << " ";
    cout << endl;
    if(flagOpt == 1)cout << "KCLUB SIZE: " << bestKClub.size() << endl;
    else cout << "KCLUB SIZE: >=" << bestKClub.size() << endl;

    if(bestKClub.size() == 0) cout << "BEST WINDOW: " << endl;
    else{
        cout << "BEST WINDOW: " << bestWindow + 1 << " (";
        for(int i = bestWindow; i < bestWindow + tauParam; i++){
            cout << "G" << i + 1;
            if(i < bestWindow + tauParam -1) cout << " ";
        }
        cout << ")" << endl;
    }
    double duration = (GetWallTime() - st);
    cout << "TIME: " << fixed << setprecision(2) << duration << " sec" << endl;

    ofstream fout;
    fout.open(methodParam + "_MW.csv", ios::in);
    if(fout.fail()){
        fout.open(methodParam + "_MW.csv", ios::out);
        fout << "instance, tau, k, obj, duration(sec)\n";
        fout.close();
    }else{
        fout.close();
    }

    fout.open(methodParam + "_MW.csv", ios::app);
    if(flagOpt == 1)
        fout << instanceName << "," << tauParam << "," << kParam << "," << bestKClub.size() << "," << fixed << setprecision(2) << duration << endl;
    else
        fout << instanceName << "," << tauParam << "," << kParam << ", >=" << bestKClub.size() << "," << fixed << setprecision(2) << duration << endl;
    fout.close();
}

void Solve(){
    double timeBegin = GetWallTime();

    vector<int> HeuSol = GetHeuristicSol();
    graph* intersectionOfPowerGraphPtr;
    bool isChanged = true;
    int edge_after_peel = 0;
    double peelTimeBegin = GetWallTime();
    //Recursive peeling--Algo 2
    while(isChanged){
        isChanged = false;
        //find total number of edges before peel
        int edge_before_peel = 0;
        for (auto g: graphCollection) {
            edge_before_peel += g->m;
        }
        intersectionOfPowerGraphPtr = GetIntersectionGraphOfPowerGraphs(); //obtain the power intersection graph J ( Gk )
        CorePeel(intersectionOfPowerGraphPtr, HeuSol.size());// core peeling
        CommunityPeel(intersectionOfPowerGraphPtr, HeuSol.size());//Community Peel
        DisconnectGraphsByComponents(graphCollection, intersectionOfPowerGraphPtr);//Cross Edge Peel

        //find total number of edges after peel
        edge_after_peel = 0;
        for (auto g: graphCollection) {
            edge_after_peel += g->m;
        }

        if (edge_before_peel != edge_after_peel)
            isChanged = true;
    }
    double peelTimeEnd = GetWallTime();
    //As the optimal solutions only exist in components whose size >= HeuSol.size()
    vector<vector<int>> componentsAfterPeeling = GetComponents(intersectionOfPowerGraphPtr, HeuSol.size());
    // get vertex set after peeling
    vector<int> verticesAfterPeeling;
    for(int i = 0; i < componentsAfterPeeling.size(); i++){
        for(int j = 0; j < componentsAfterPeeling[i].size(); j++){
            verticesAfterPeeling.push_back(componentsAfterPeeling[i][j]);
        }
    }
    //sort
    sort(verticesAfterPeeling.begin(),verticesAfterPeeling.end());
    // get intersection of power graph after peeling
    graph* intersectionOfPowerGraphPtrAfterPeeling = intersectionOfPowerGraphPtr->CreateSubgraph(&verticesAfterPeeling);
    vector<graph*> graphCollectionAfterPeeling;

    //get number of vertices and edges after peeling
    num_vertex_after_peel = verticesAfterPeeling.size();
    num_edges_collection_after_peel = 0;
    graph *g;
    for(auto i : graphCollection){
        g = i->CreateSubgraph(&verticesAfterPeeling);
        num_edges_collection_after_peel += g->m;
        graphCollectionAfterPeeling.push_back(g);
    }
    //get components of power intersection graph
    componentsAfterPeeling = GetComponents(intersectionOfPowerGraphPtrAfterPeeling, 0);

    double timeMiddle = GetWallTime();
    //got optimal solution through Heuristic
    if(componentsAfterPeeling.size() == 0){
        cout << "\nSUMMARY" << endl;
        cout << "PARAMETER P : " << pParam << endl;
        cout << "PARAMETER K : " << kParam << endl;
        cout << "INSTANCE : " << instanceName << " " << firstGraph << "-" << lastGraph << endl;
        cout << "METHOD : " << methodParam << endl;
        cout << "number of nodes before peeling: "<<num_vertex_before_peel << endl;
        cout << "number of nodes after peeling: "<<num_vertex_after_peel << endl;
        cout << "number of edges before peeling: "<<num_edges_collection_before_peel << endl;
        cout << "number of edges after peeling: "<<num_edges_collection_after_peel << endl;
        cout << "OBJ VALUE : " << HeuSol.size() << endl;
        cout<<"Optimal solution by Heuristic is: **********************"<<endl;
        for (int i = 0; i < HeuSol.size(); i++) {
            if (i<HeuSol.size()-1){
                cout<<HeuSol[i] +1 <<"->";
            } else
                cout<<HeuSol[i] +1<<endl;
        }
        cout << "OBJ BOUND : " << "-" << endl;
        cout << "DURATION : " << fixed << setprecision(2) << timeMiddle - timeBegin << endl;
        cout << "NUM OF CALLBACKS : " << "-" << endl;
        cout << "CALLBACK TIME : " << "-" << endl;
        cout << "NUM OF LAZY CONSTRAINTS ADDED : " << "-" << endl;
        cout << "NUM OF STRENGTHENED LAZY CONSTRAINTS ADDED : " << "-" << endl;
        cout << "AVE NUM OF NEGATIVE TERMS : " << "-" << endl;
        cout << "AVE NUM OF NEGATIVE TERMS STRENGTHENED : " << "-" << endl;
        cout << "NUM OF BB NODES : " << "-" << endl;

        ofstream fout;
        fout.open(methodParam + ".csv", ios::in);
        if(fout.fail()){
            fout.open(methodParam + ".csv", ios::out);
            fout << "method,p, k, instance,#nodes before peel,#nodes after peel,#edges before peel, #edge after peel, Heuristic value, obj, obj bound, duration, num callbacks, callback time, num lazy cuts, num strengthened lazy cuts, ave num of negative terms, ave num of negative terms strengthened, num BB nodes\n";
            fout.close();
        }else{
            fout.close();
        }

        fout.open(methodParam + ".csv", ios::app);
        fout << methodParam << ",";
        fout << pParam << ",";
        fout << kParam << ",";
        fout << instanceName << " " << firstGraph << "-" << lastGraph << ",";

        fout << num_vertex_before_peel << ",";
        fout << num_vertex_after_peel << ",";
        fout << num_edges_collection_before_peel << ",";
        fout << num_edges_collection_after_peel << ",";
        fout << HeuSol.size() << ","; //heuristic size
        fout << HeuSol.size() << ","; //objective value (because here heuristic is optimal solution)
        fout << "-" << ",";
        fout << fixed << setprecision(2) << timeMiddle - timeBegin << ",";
        fout << "-" << ",";
        fout << "-" << ",";
        fout << "-" << ",";
        fout <<  "-" << ",";
        fout << "-" << ",";
        fout << "-" << ",";
        fout << "-" << endl;
        fout.close();
        cout << "\nFINISH WRITING OUT RESULTS" << endl;
        return;
    }
    //MIP model
    try{

        for(auto graphPtr : graphCollectionAfterPeeling)
            graphPtr->FindKNeighbors(kParam);

        GRBEnv* ENV = 0;
        GRBVar* X = 0;
        GRBVar* Y = 0;

        ENV = new GRBEnv();
        GRBModel MODEL = GRBModel(*ENV);
        int nParam = graphCollectionAfterPeeling[0]->n;

        X = MODEL.addVars(nParam, GRB_BINARY);
        Y = MODEL.addVars(componentsAfterPeeling.size(), GRB_BINARY);
        MODEL.update();

        for(int i = 0; i < nParam; i++)
            X[i].set(GRB_DoubleAttr_Obj, 1.0);
        MODEL.update();

        GRBLinExpr sumY;
        for(int i = 0; i < componentsAfterPeeling.size(); i++){
            sumY += Y[i];
            for(int j = 0; j < componentsAfterPeeling[i].size(); j++){
                MODEL.addConstr(X[componentsAfterPeeling[i][j]] <= Y[i]);
            }
        }
        MODEL.addConstr(sumY == 1);

        for(int i = 0; i < componentsAfterPeeling.size(); i++){
            for(int u = 0; u < componentsAfterPeeling[i].size(); u++){
                int uVertex = componentsAfterPeeling[i][u];
                for(int v = u + 1; v < componentsAfterPeeling[i].size(); v++){
                    int vVertex = componentsAfterPeeling[i][v];
                    if(intersectionOfPowerGraphPtrAfterPeeling->IsAdjacent(uVertex,vVertex) == false)
                        MODEL.addConstr(X[uVertex] + X[vVertex] <= 1);
                }
            }
        }
        MODEL.update();
        /*
        MODEL.set(GRB_IntParam_MIPFocus, 0); //0 default value
        MODEL.set(GRB_IntParam_Threads, 0); //0 default value
        MODEL.set(GRB_IntParam_Method, -1);//Algorithm used to solve continuous models, -1=automatic,
        MODEL.set(GRB_IntParam_NodeMethod, 1);//Algorithm used for MIP node relaxations (except for the initial root node relaxation, see Method). Options are: -1=automatic, 0=primal simplex, 1=dual simplex, and 2=barrier.
        */
         //set global cut aggressiveness; over-ridden by individual cut settings
        MODEL.set(GRB_IntParam_Cuts, 0);
        //0=no cuts;1=moderate;2=aggressive;3=very aggressive;-1=default
        MODEL.set(GRB_DoubleParam_TimeLimit, 7200);
        //Set termination gap limit; as needed; default is 1e-4
        //MODEL.set(GRB_DoubleParam_MIPGap, 1e-4);
        MODEL.set(GRB_StringParam_LogFile, "");
        MODEL.set(GRB_IntParam_OutputFlag, 1); // turn on log printing
        MODEL.set(GRB_IntAttr_ModelSense, -1);//objective maximization
        MODEL.set(GRB_IntParam_LazyConstraints,1);
        int numStrengthenedLazyCut = 0, numLazyCut = 0, numCallback = 0, independentSetInequalityCount = 0, independentSetInequalityCountViolateEpsilon = 0;
        double callbackTime = 0, aveNumTermsNeg = 0, aveNumTermsNegStrengthened = 0;

        vector<GRBModel> SEP_MODELS; //if Method does not use PPCF_IS, they are empty set; this is to have a general callback function
        vector<GRBVar*> Z(pParam);
        vector<GRBVar*> W(pParam);
        if(methodParam == "PPCF_IS"){
            // use PPCF + Independent set inequality
            // Separation formulation begin
            int nParam_Sep = graphCollectionAfterPeeling[0]->n;
            vector<GRBEnv*> ENVS_SEP(pParam);
            for(int i = 0; i < pParam; i++){
                ENVS_SEP[i] = 0;
                Z[i] = 0;
                W[i] = 0;
                ENVS_SEP[i] = new GRBEnv();

                GRBModel MODEL_SEP = GRBModel(*(ENVS_SEP[i]));
                MODEL_SEP.set(GRB_IntAttr_ModelSense, -1);

                MODEL_SEP.set(GRB_DoubleParam_TimeLimit, 30);
                //Set Gurobi screen display flag: 0=switch off; 1=default
                MODEL_SEP.getEnv().set(GRB_IntParam_OutputFlag,0);
                Z[i] = MODEL_SEP.addVars(nParam_Sep, GRB_BINARY);
                W[i] = MODEL_SEP.addVars(nParam_Sep, GRB_CONTINUOUS);
                MODEL_SEP.update();

                for(int u = 0; u < nParam_Sep; u++){
                    W[i][u].set(GRB_DoubleAttr_LB, 0.0);
                    W[i][u].set(GRB_DoubleAttr_Obj, 0.0);
                    Z[i][u].set(GRB_DoubleAttr_Obj, 0.0);

                    stringstream ss;
                    ss << u;
                    W[i][u].set(GRB_StringAttr_VarName, "W" + ss.str());
                    Z[i][u].set(GRB_StringAttr_VarName, "Z" + ss.str());

                    MODEL_SEP.update();
                    vector<int> temp = FindCommon(graphCollectionAfterPeeling, i, u);
                    MODEL_SEP.addConstr(W[i][u] <= (temp.size())*(1 - Z[i][u]));

                    GRBLinExpr RHS = 0;
                    for(int w = 0; w < temp.size(); w++){
                        RHS += Z[i][temp[w]];
                    }
                    RHS = RHS - 1 - (temp.size())*Z[i][u];
                    MODEL_SEP.addConstr(W[i][u] >= RHS);

                    for(int v = u + 1; v < nParam_Sep; v++){
                        if(graphCollectionAfterPeeling[i]->IsAdjacent(u, v)){
                            MODEL_SEP.addConstr(Z[i][u] + Z[i][v] <= 1);
                        }
                    }
                }
                MODEL_SEP.update();
                SEP_MODELS.push_back(MODEL_SEP);
            }
            //Separation formulation end
        }//end if
        pGraph_callback myCallback = pGraph_callback(methodParam, X, graphCollectionAfterPeeling, kParam, numStrengthenedLazyCut, numLazyCut, numCallback, callbackTime, aveNumTermsNeg, aveNumTermsNegStrengthened, SEP_MODELS, Z, W, epsilonParam, independentSetInequalityCount, independentSetInequalityCountViolateEpsilon);
        MODEL.setCallback(&myCallback);
        MODEL.optimize();
        double timeEnd = GetWallTime();

        //validate solution
        //MODEL.write("model.lp");
        if(MODEL.get(GRB_IntAttr_SolCount)){
            vector<int>* solutionPtr = new vector<int>();
            cout<<"The optimal solution (with original node ID) is:---"<<endl;
            for(int i = 0; i < nParam; i++){
                if(X[i].get(GRB_DoubleAttr_X) > 0.5){
                    solutionPtr->push_back(intersectionOfPowerGraphPtrAfterPeeling->vertices[i]->name);
                    cout<<intersectionOfPowerGraphPtrAfterPeeling->vertices[i]->name + 1 << "->";
                    //cout << i << " " << intersectionOfPowerGraphPtrAfterPeeling->vertices[i]->name << " " << X[i].get(GRB_StringAttr_VarName) << endl;
                }
            }
            cout<<endl;

            if(ValidateSolution(solutionPtr) == false){
                delete solutionPtr;
                return;
            }
            delete solutionPtr;
        }else if(HeuSol.size() > 0){
            if(ValidateSolution(&HeuSol) == false){
                return;
            }
        }
        else{
            cout << "NO SOLUTION FOUND SO FAR" << endl;
            return;
        }

        //cout and write out results
        double objValue = MODEL.get(GRB_DoubleAttr_ObjVal);
        double objBound = MODEL.get(GRB_DoubleAttr_ObjBound);
        double MIP_Gap = 100*(objBound-objValue)/objValue;
        int numBBNode = MODEL.get(GRB_DoubleAttr_NodeCount);
        double PeelTime = peelTimeEnd - peelTimeBegin;
        double duration = timeEnd - timeBegin;
        double grbSolveTime =  MODEL.get(GRB_DoubleAttr_Runtime);

        delete ENV;

        //cout << green << "\nSUMMARY" << def << endl;
        cout << "\nSUMMARY" << endl;
        cout << "PARAMETER P : " << pParam << endl;
        cout << "PARAMETER K : " << kParam << endl;
        cout << "INSTANCE : " << instanceName << " " << firstGraph << "-" << lastGraph << endl;
        cout << "METHOD : " << methodParam << endl;
        cout << "number of nodes before peeling: "<<num_vertex_before_peel << endl;
        cout << "number of nodes after peeling: "<<num_vertex_after_peel << endl;
        cout << "number of edges before peeling: "<<num_edges_collection_before_peel << endl;
        cout << "number of edges after peeling: "<<num_edges_collection_after_peel << endl;
        cout << "OBJ VALUE : " << objValue << endl;
        cout << "OBJ BOUND : " << objBound << endl;
        cout << "MIP Gap (%) : " << MIP_Gap << endl;
        cout << "Peel Time : " << PeelTime << endl;
        cout << "Gurobi solve Time : " << grbSolveTime<< endl;
        cout << "DURATION : " << fixed << setprecision(2) << duration << endl;
        cout << "NUM OF CALLBACKS : " << numCallback << endl;
        cout << "CALLBACK TIME : " << fixed << setprecision(2) << callbackTime << endl;
        cout << "NUM OF LAZY CONSTRAINTS ADDED : " << numLazyCut << endl;
        cout << "NUM OF STRENGTHENED LAZY CONSTRAINTS ADDED : " << numStrengthenedLazyCut << endl;
        cout << "AVE NUM OF NEGATIVE TERMS : " << aveNumTermsNeg << endl;
        cout << "AVE NUM OF NEGATIVE TERMS STRENGTHENED : " << aveNumTermsNegStrengthened << endl;
        cout << "NUM OF BB NODES : " << numBBNode << endl;

        ofstream fout;
        fout.open(methodParam + ".csv", ios::in);
        if(fout.fail()){
            fout.open(methodParam + ".csv", ios::out);
            fout << "method,p, k, instance,#nodes before peel,#nodes after peel,#edges before peel, #edge after peel,Heuristic value, obj, obj bound,MIP Gap (%), Peel Time, GRB solve Time, duration, num callbacks, callback time, num lazy cuts, num strengthened lazy cuts, ave num of negative terms, ave num of negative terms strengthened, num BB nodes, num ind set cuts violated, num ind set cuts violated epsilon\n";
            fout.close();
        }else{
            fout.close();
        }

        fout.open(methodParam + ".csv", ios::app);
        fout << methodParam << ",";
        fout << pParam << ",";
        fout << kParam << ",";
        fout << instanceName << " " << firstGraph << "-" << lastGraph << ",";
        fout << num_vertex_before_peel << ",";
        fout << num_vertex_after_peel << ",";
        fout << num_edges_collection_before_peel << ",";
        fout << num_edges_collection_after_peel << ",";
        fout << HeuSol.size() << ",";
        fout << objValue << ",";
        fout << objBound << ",";
        fout << MIP_Gap << ",";
        fout << PeelTime << ",";
        fout << grbSolveTime << ",";
        fout << fixed << setprecision(2) << duration << ",";
        fout << numCallback << ",";
        fout << fixed << setprecision(2) << callbackTime << ",";
        fout << numLazyCut << ",";
        fout <<  numStrengthenedLazyCut << ",";
        fout << aveNumTermsNeg << ",";
        fout << aveNumTermsNegStrengthened << ",";
        fout << numBBNode << ",";
        fout << independentSetInequalityCount << ",";
        fout << independentSetInequalityCountViolateEpsilon << endl;
        fout.close();
        //cout << green << "\nFINISH WRITING OUT RESULTS" << def << endl;
        cout << "\nFINISH WRITING OUT RESULTS" << endl;
    }catch (GRBException e) {
        cout << "ERROR CODE = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "EXCEPTION DURING OPTIMIZATION" << endl;
    }

    for(auto i : graphCollectionAfterPeeling) delete i;
    delete intersectionOfPowerGraphPtr;
    delete intersectionOfPowerGraphPtrAfterPeeling;
}


//Cross-edge peeling

void DisconnectGraphsByComponents(vector<graph*> graphCollection, graph* intersectionOfPowerGraphPtr){

    graphCollection.push_back(intersectionOfPowerGraphPtr);

    int numEdgesDeleted;
    do{
        numEdgesDeleted = 0;
        vector<vector<int>> componentList;
        for (auto g: graphCollection) {
            componentList = GetComponents(g, 0);
            //Examine each pair of connected components within each graph, verifying whether there exist nodes u and v  in distinct components, yet forming an edge in some graph.
            int size_compList = int(componentList.size());
            //if the number of  components are equal to number of vertices, it means that all vertices are isolated; then we can reset every graph
            if (size_compList == g->n){
                //reset each graph to empty
                for(auto g2: graphCollection){
                    g2->initialize();
                }
                return;
            }
            for (int i = 0; i < size_compList; i++) {
                for (int j = i+1; j < size_compList; j++) {
                    DeleteEdgesBasedOnComponents(graphCollection, componentList[i], componentList[j], numEdgesDeleted);
                }
            }
        }
    }while(numEdgesDeleted > 0);

    // update max degree
    for(auto i : graphCollection){
        int maxDeg = 0;
        int maxDegreeVertexIndex = -1;
        for(int j = 0; j < i->n; j++){
            if(i->vertices[j]->degree > maxDeg){
                maxDeg = (int)i->vertices[j]->degree;
                maxDegreeVertexIndex = j;
            }
        }
        i->maxDegree = maxDeg;
        i->maxDegreeVertexIndex = maxDegreeVertexIndex;
    }

}



void DeleteEdgesBasedOnComponents(vector<graph*> graphCollection, vector<int>& compOne, vector<int>& compTwo, int& numEdgesDeleted){
    //CompOne and CompOne are in the same graph, so they do have any overlap
    if (compOne.empty() && compTwo.empty()) return;
    for (auto i:compOne) {
        for(auto j: compTwo){
            for(auto graphPtr : graphCollection){
                if(graphPtr->IsAdjacent(i, j)){
                    graphPtr->vertices[i]->neighbors.erase(find(graphPtr->vertices[i]->neighbors.begin(), graphPtr->vertices[i]->neighbors.end(), j));
                    graphPtr->vertices[j]->neighbors.erase(find(graphPtr->vertices[j]->neighbors.begin(), graphPtr->vertices[j]->neighbors.end(), i));
                    graphPtr->vertices[i]->degree--;
                    graphPtr->vertices[j]->degree--;
                    graphPtr->m--;
                    numEdgesDeleted++;
                }
            }
        }
    }
}


//core peel power graph + graphs collections
void CorePeel(graph* graphP, int peel){
    if(graphP->maxDegree > 0){
        queue<int> Q;
        vector<int> F(graphP->n, 0);
        for(int i = 0; i < graphP->n; i++){
            if(graphP->vertices[i]->degree < peel){
                Q.push(i);
                F[i] = 1;
            }
        }

        while(Q.size()){
            int pVertex = Q.front();
            Q.pop();
            graphP->vertices[pVertex]->degree = 0;
            for(int i = 0; i < graphP->vertices[pVertex]->neighbors.size(); i++){
                int pNeighbor = graphP->vertices[pVertex]->neighbors[i];
                if(graphP->vertices[pNeighbor]->degree > 0){
                    graphP->vertices[pNeighbor]->degree--;
                    if(graphP->vertices[pNeighbor]->degree > 0){
                        if(graphP->vertices[pNeighbor]->degree < peel){
                            if(F[pNeighbor] == 0){
                                Q.push(pNeighbor);
                                F[pNeighbor] = 1;
                            }
                        }
                    }
                }
            }
        }

        vector<int> neighbors;
        int maxDeg = 0;
        int maxDegNode = -1;
        int edgeCount = 0;
        for(int i = 0; i < graphP->n; i++){
            neighbors.clear();
            if(graphP->vertices[i]->degree > 0){
                for(int j = 0; j < graphP->vertices[i]->neighbors.size(); j++){
                    int pNeighbor = graphP->vertices[i]->neighbors[j];
                    if(graphP->vertices[pNeighbor]->degree > 0){
                        neighbors.push_back(pNeighbor);
                    }
                }
                graphP->vertices[i]->neighbors = neighbors;
                edgeCount += (int)neighbors.size();
                if(graphP->vertices[i]->degree > maxDeg){
                    maxDeg = graphP->vertices[i]->degree;
                    maxDegNode = i;
                }

            }else{
                graphP->vertices[i]->neighbors.clear();
            }
        }
        graphP->m = edgeCount/2;
        if(maxDeg > 0){
            graphP->maxDegree = maxDeg;
            graphP->maxDegreeVertexIndex = maxDegNode;
        }else{
            graphP->maxDegree = 0;
            graphP->maxDegreeVertexIndex = 0;
        }

        //update Graph collection: if a vertex's degree =0 (equivalent to the deletion of such vertex), then such vertex should be removed from every graph in the collection

        vector<bool> isVertexRemoved(graphP->n, false); //find vertices to be removed; set true if yes
        for(int i = 0; i < graphP->n; i++){
            if (graphP->vertices[i]->degree ==0)
                isVertexRemoved[i] = true;
        }

        for(auto g : graphCollection){
            //removed vertices
            g->GraphRemoveNodes(isVertexRemoved);
        }
    }
}

//Community Peel powergraph + graph collection

void CommunityPeel(graph* graphP, int peel){
    if(graphP->maxDegree > 0){
        queue<vector<int>> Q;
        //initialization
        for(int i = 0; i < graphP->n; i++){
            graphP->vertices[i]->adjacency.resize(graphP->n, 0);
            graphP->vertices[i]->F.resize(graphP->n, 0); // used as flag to indicate whether an edge is put into Q or not for peeling
            graphP->vertices[i]->numCommonNeighbors.resize(graphP->n, 0);
            graphP->vertices[i]->commonNeighbors.resize(graphP->n);
        }
        //1. find pairs of edges to be peeling
        for(int i = 0; i < graphP->n; i++){
            for(int j = 0; j < graphP->vertices[i]->neighbors.size(); j++){
                int pNeighbor = graphP->vertices[i]->neighbors[j];
                if(graphP->vertices[i]->adjacency[pNeighbor] == 0){
                    graphP->vertices[i]->adjacency[pNeighbor] = 1;
                    graphP->vertices[pNeighbor]->adjacency[i] = 1;

                    graphP->vertices[i]->commonNeighbors[pNeighbor] = graphP->FindCommonNeighbors(i, pNeighbor);
                    graphP->vertices[i]->numCommonNeighbors[pNeighbor] = (int)graphP->vertices[i]->commonNeighbors[pNeighbor].size();
                    graphP->vertices[pNeighbor]->commonNeighbors[pNeighbor] = graphP->vertices[i]->commonNeighbors[pNeighbor];
                    graphP->vertices[pNeighbor]->numCommonNeighbors[i] = graphP->vertices[i]->numCommonNeighbors[pNeighbor];

                    if(graphP->vertices[i]->numCommonNeighbors[pNeighbor] < (peel - 1)){
                        if(graphP->vertices[i]->F[pNeighbor] == 0){
                            vector<int> tempVec{i, pNeighbor};
                            Q.push(tempVec);
                            graphP->vertices[i]->F[pNeighbor] = 1;
                            graphP->vertices[pNeighbor]->F[i] = 1;
                        }
                    }
                }
            }
        }
        //recursive edge peeling
        while(Q.size()){
            vector<int> pEdge = Q.front();
            Q.pop();
            int u = pEdge[0], v = pEdge[1];
            int i, j;
            if(u < v){
                i = u;
                j = v;
            }else{
                i = v;
                j = u;
            }

            graphP->vertices[i]->adjacency[j] = 0;
            graphP->vertices[j]->adjacency[i] = 0;
            for(int l = 0; l < graphP->vertices[i]->commonNeighbors[j].size(); l++){
                int pCommonNeighbor = graphP->vertices[i]->commonNeighbors[j][l];

                if(graphP->vertices[i]->adjacency[pCommonNeighbor] > 0 && graphP->vertices[j]->adjacency[pCommonNeighbor] > 0){
                    graphP->vertices[i]->numCommonNeighbors[pCommonNeighbor]--;
                    graphP->vertices[pCommonNeighbor]->numCommonNeighbors[i]--;
                    if(graphP->vertices[i]->numCommonNeighbors[pCommonNeighbor] < (peel - 1)){
                        if(graphP->vertices[i]->F[pCommonNeighbor] == 0){
                            vector<int> tempVec{i, pCommonNeighbor};
                            Q.push(tempVec);
                            graphP->vertices[i]->F[pCommonNeighbor] = 1;
                            graphP->vertices[pCommonNeighbor]->F[i] = 1;
                        }
                    }

                    graphP->vertices[j]->numCommonNeighbors[pCommonNeighbor]--;
                    graphP->vertices[pCommonNeighbor]->numCommonNeighbors[j]--;
                    if(graphP->vertices[j]->numCommonNeighbors[pCommonNeighbor] < (peel - 1)){
                        if(graphP->vertices[j]->F[pCommonNeighbor] == 0){
                            vector<int> tempVec{j, pCommonNeighbor};
                            Q.push(tempVec);
                            graphP->vertices[j]->F[pCommonNeighbor] = 1;
                            graphP->vertices[pCommonNeighbor]->F[j] = 1;
                        }
                    }
                }

            }
        }
        //update graphP
        vector<int> neighbors;
        int maxDeg = 0;
        int maxDegNode = -1;
        int edgeCount = 0;
        for(int i = 0; i < graphP->n; i++){
            neighbors.clear();
            for(int j = 0; j < graphP->vertices[i]->neighbors.size(); j++){
                int pNeighbor = graphP->vertices[i]->neighbors[j];
                if(graphP->vertices[i]->adjacency[pNeighbor] == 1){
                    neighbors.push_back(pNeighbor);
                }
            }
            graphP->vertices[i]->neighbors = neighbors;
            graphP->vertices[i]->degree = (int)neighbors.size();
            edgeCount += graphP->vertices[i]->degree;
            if(maxDeg < graphP->vertices[i]->degree){
                maxDeg = graphP->vertices[i]->degree;
                maxDegNode = i;
            }
        }
        graphP->m = edgeCount/2;
        if(maxDeg > 0){
            graphP->maxDegree = maxDeg;
            graphP->maxDegreeVertexIndex = maxDegNode;
        }else{
            graphP->maxDegree = 0;
            graphP->maxDegreeVertexIndex = 0;
        }

        //update graph collection by removing edges
        //1. find pairs of edges {i,j} that were deleted in power graph--graphP
        vector<vector<int>> edge_removal_list;
        for(int i = 0; i < graphP->n; i++){
            for(int j = i + 1; j < graphP->n; j++){
                //2. remove such edge {i,j} from every graph in the collection if exist
                if (graphP->vertices[i]->F[j] == 1){
                    for(auto g : graphCollection){
                        g->DeleteEdge(i,j); //remove edge
                    }
                }
            }
        }

        //variables clear
        for(int i = 0; i < graphP->n; i++){
            vector<int>().swap(graphP->vertices[i]->adjacency);
            vector<int>().swap(graphP->vertices[i]->F);
            vector<int>().swap(graphP->vertices[i]->numCommonNeighbors);
            vector<vector<int>>().swap(graphP->vertices[i]->commonNeighbors);
        }
    }
}


// get components of a graph
vector<vector<int>> GetComponents(graph* graphP, int peel){
    //peel is the minimum component size
    vector<vector<int>> components;
    vector<int> F(graphP->n, 0);
    vector<int> component;
    for(int i = 0; i < graphP->n; i++){
        if(F[i] == 0){
            component.clear();
            F[i] = 1;
            component = graphP->BFS(i);
            component.push_back(i);
            sort(component.begin(), component.end());
            for(int j = 0; j < component.size(); j++){
                int pNode = component[j];
                F[pNode] = 1;
            }
            if((int)component.size() >= peel){
                components.push_back(component);
            }
        }
    }

    vector<int>().swap(F);
    vector<int>().swap(component);
    sort(components.begin(), components.end(), cmp);// sort components in descending order of component size
    return components;
}


// get intersection of power graphs
graph* GetIntersectionGraphOfPowerGraphs(){

    for (int i = 0; i < graphCollection.size(); i++) {
        graphCollection[i]->FindKNeighbors(kParam);
    }

    graph* graphPtr = new graph();
    graphPtr->n = graphCollection[0]->n;

    for(int i = 0; i < graphPtr->n; i++){
        vertex* vertexPtr = new vertex();
        vertexPtr->ind = i;
        graphPtr->vertices.push_back(vertexPtr);
    }

    int edgeCount = 0;
    for(int i = 0; i < graphPtr->n; i++){
        for(int j = 0; j < graphCollection[0]->vertices[i]->kNeighbors.size(); j++){
            int pNeighbor = graphCollection[0]->vertices[i]->kNeighbors[j];
            if(pNeighbor > i){
                int graphCount = 1;
                for(int t =  1; t < graphCollection.size(); t++){
                    if(graphCollection[t]->IsKAdjacent(i, pNeighbor)){
                        graphCount++;
                    }else{
                        graphCount = 0;
                        break;
                    }
                }
                if(graphCount > 0){
                    graphPtr->vertices[i]->neighbors.push_back(pNeighbor);
                    graphPtr->vertices[pNeighbor]->neighbors.push_back(i);
                    edgeCount++;
                }
            }
        }
    }
    graphPtr->m = edgeCount;

    int maxDeg = 0;
    int maxDegNode = -1;
    for(int i = 0; i < graphPtr->n; i++){
        graphPtr->vertices[i]->degree = graphPtr->vertices[i]->neighbors.size();
        if(graphPtr->vertices[i]->degree > maxDeg){
            maxDeg = graphPtr->vertices[i]->degree;
            maxDegNode = i;
        }
    }
    graphPtr->maxDegree = maxDeg;
    graphPtr->maxDegreeVertexIndex = maxDegNode;
    return graphPtr;
}


vector<vector<int>> GetComponents(graph* graphPtr){
    vector<vector<int>> components;
    vector<int> F(graphPtr->n, 0);
    vector<int> component;
    for(int i = 0; i < graphPtr->n; i++){
        if(F[i] == 0){
            component.clear();
            F[i] = 1;
            component = graphPtr->BFS(i);
            component.push_back(i);
            sort(component.begin(), component.end());
            for(int j = 0; j < component.size(); j++){
                int pVertex = component[j];
                F[pVertex] = 1;
            }

            components.push_back(component);

        }

        sort(components.begin(), components.end(), cmp);// sort components in descending order of component size
    }

    return components;
}





vector<int>* MINIMALIZE(graph* graphPtr, int u, int v, int kParam, vector<int>* trivialSeparator){

    vector<bool>* minimalSeparatorBool = new vector<bool>(graphPtr->n, false);
    graphPtr->CalculateDistanceFrom(u);
    graphPtr->CalculateDistanceFrom(v);

    for(int i = 0; i < trivialSeparator->size(); i++){
        int iVertex = (*trivialSeparator)[i];
        if(graphPtr->vertices[u]->distanceTo[iVertex] + graphPtr->vertices[v]->distanceTo[iVertex] <= kParam)
            (*minimalSeparatorBool)[iVertex] = true;
    }

    vector<bool>* minimalSeparatorBoolC = new vector<bool>(graphPtr->n, true);
    for(int i = 0; i < graphPtr->n; i++){
        if((*minimalSeparatorBool)[i])
            (*minimalSeparatorBoolC)[i] = false;
    }

    for(int i = 0; i < graphPtr->n; i++){

        if((*minimalSeparatorBool)[i] == false) continue;
        if(graphPtr->vertices[u]->distanceTo[i] == 1 && graphPtr->vertices[v]->distanceTo[i] == 1) continue;
        (*minimalSeparatorBoolC)[i] = true;

        graphPtr->CalculateDistanceFromTo(i, minimalSeparatorBoolC);

        if(graphPtr->vertices[i]->distanceTo[u] + graphPtr->vertices[i]->distanceTo[v] <= kParam)
            (*minimalSeparatorBoolC)[i] = false;
        else
            (*minimalSeparatorBool)[i] = false;
    }

    vector<int>* minimalSeparatorPtr = new vector<int>();
    for(int i = 0; i < graphPtr->n; i++){
        if((*minimalSeparatorBool)[i])
            minimalSeparatorPtr->push_back(i);
    }

    delete minimalSeparatorBool;
    delete minimalSeparatorBoolC;

    return minimalSeparatorPtr;
}


vector<int>* MINIMALIZE(graph* graphPtr, int u, int v, int kParam, vector<int>* toBeDeleted, vector<int>* trivialSeparator){

    vector<int>* minimalSeparatorPtr = new vector<int>();
    minimalSeparatorPtr->clear();

    vector<bool>* toBeKeptBool = new vector<bool>(graphPtr->n, true);
    for(auto i : (*toBeDeleted))
        (*toBeKeptBool)[i] = false;

    graphPtr->CalculateDistanceFromTo(u, toBeKeptBool);
    graphPtr->CalculateDistanceFromTo(v, toBeKeptBool);

    if(graphPtr->vertices[u]->distanceTo[v] > kParam){
        delete toBeKeptBool;
        return minimalSeparatorPtr;
    }

    vector<bool>* minimalSeparatorBool = new vector<bool>(graphPtr->n, false);
    for(auto i : (*trivialSeparator)){
        if(graphPtr->vertices[u]->distanceTo[i] + graphPtr->vertices[v]->distanceTo[i] <= kParam)
            (*minimalSeparatorBool)[i] = true;
    }

    for(int i = 0; i < graphPtr->n; i++){
        if((*minimalSeparatorBool)[i])
            (*toBeKeptBool)[i] = false;
    }

    for(int i = 0; i < graphPtr->n; i++){

        if((*minimalSeparatorBool)[i] == false) continue;
        if(graphPtr->vertices[u]->distanceTo[i] == 1 && graphPtr->vertices[v]->distanceTo[i] == 1) continue;
        (*toBeKeptBool)[i] = true;

        graphPtr->CalculateDistanceFromTo(i, toBeKeptBool);

        if(graphPtr->vertices[i]->distanceTo[u] + graphPtr->vertices[i]->distanceTo[v] <= kParam)
            (*toBeKeptBool)[i] = false;
        else
            (*minimalSeparatorBool)[i] = false;
    }


    for(int i = 0; i < graphPtr->n; i++){
        if((*minimalSeparatorBool)[i])
            minimalSeparatorPtr->push_back(i);
    }

    delete toBeKeptBool;
    delete minimalSeparatorBool;

    return minimalSeparatorPtr;
}



bool ValidateSolution(vector<int>* solutionPtr){

    graph* subgraphPtr;

    for(auto graphPtr : graphCollection){
        subgraphPtr = graphPtr->CreateSubgraph(solutionPtr);
        for(int i = 0; i < subgraphPtr->n; i++){
            subgraphPtr->CalculateDistanceFrom(i);
            for(int j = i + 1; j < subgraphPtr->n; j++){
                if(subgraphPtr->vertices[i]->distanceTo[j] > kParam){
                    //cout << green << "\nINVALID SOLUTION DETECTED" << def << endl;
                    cout << "\nINVALID SOLUTION DETECTED" << endl;
                    cout << "INVALID VERTEX PAIR : " << subgraphPtr->vertices[i]->name << " " << subgraphPtr->vertices[j]->name << endl;
                    delete subgraphPtr;
                    return false;
                }
            }
        }
        delete subgraphPtr;
    }

    //cout << green << "\nSOLUTION VALIDATED" << def << endl;
    cout << "\nSOLUTION VALIDATED" << endl;
    return true;
}


void CleanUp(){
    for(auto graphPtr : graphCollection)
        delete graphPtr;
    graphCollection.clear();
    //cout << green << "\nMEMORY DEALLOCATED FOR GRAPH COLLECTION" << def << endl;
    cout << "\nMEMORY DEALLOCATED FOR GRAPH COLLECTION" << endl;
}

void CleanUpKClubSig(){
    for(auto graphPtr : graphSequence)
        delete graphPtr;
    graphSequence.clear();
    cout << "\nMEMORY DEALLOCATED FOR GRAPH SEQUENCE" << endl;
}


void CalculateRatio(){

    vector<double>* theRatios = new vector<double>();
    int nParam = graphCollection[0]->n;

    for(auto graphPtr : graphCollection)
        graphPtr->FindKNeighbors(kParam);

    for(int i = 0; i < nParam; i++){
        for(int j = i + 1; j < nParam; j++){
            //cout << IsKAdjacentOverCollection(i,j) << endl;
            if((!IsAdjacentOverCollection(i,j)) && IsKAdjacentOverCollection(i,j)){

                vector<vector<int>*>* common = new vector<vector<int>*>();

                for(auto graphPtr : graphCollection)
                    common->push_back(graphPtr->FindCommonKNeighbors(i, j));

                vector<int>* commonCommon = TakeIntersection(common);

                double tempSum = 0;
                for(auto k : (*common)) tempSum += k->size();
                double tempRatio = commonCommon->size() * graphCollection.size() / tempSum;

                theRatios->push_back(tempRatio);

                delete commonCommon;
                for(auto k : (*common)) delete k;
                delete common;
            }
        }
    }

    double theRatio = 0;

    if((*theRatios).size() > 0){
        for(auto tempRatio : (*theRatios)) {theRatio += tempRatio;} //cout << tempRatio << " " << theRatio << endl;}
        theRatio = theRatio / (*theRatios).size();
    }else{
        theRatio = -1;
    }

    delete theRatios;

    fstream fout;
    fout.open("RATIO.csv", std::fstream::in | std::fstream::out | std::fstream::app);
    if(!fout){
        fout.open("RATIO.csv", fstream::in | fstream::out | fstream::trunc);
        fout << "\n";
        fout.close();
    }else{
        fout.close();
    }

    fout.open("RATIO.csv", std::fstream::in | std::fstream::out | std::fstream::app);
    fout << pParam << ",";
    fout << kParam << ",";
    fout << instanceName << ",";
    fout << theRatio << endl;
    fout.close();

    cout << "RATIO CALCULATED" << endl;
}


bool IsAdjacentOverCollection(int u, int v){

    bool answer = true;

    for(auto graphPtr : graphCollection){
        if(graphPtr->IsAdjacent(u, v)) continue;
        else{
            answer = false;
            break;
        }
    }
    return answer;
}


bool IsKAdjacentOverCollection(int u, int v){

    bool answer = true;

    for(auto graphPtr : graphCollection){
        if(graphPtr->IsKAdjacent(u, v)) continue;
        else{
            answer = false;
            break;
        }
    }
    return answer;
}


vector<int>* TakeIntersection(vector<vector<int>*>* common){

vector<int>* theIntersection = new vector<int>();
theIntersection->clear();

vector<int>* positions = new vector<int>((*common).size(), 0);

bool flag = true;
for(int i = 0; i < positions->size(); i++){
if((*positions)[i] < (*common)[i]->size()) continue;
else{
flag = false;
break;
}
}

while(flag){

set<int>* temp = new set<int>();

for(int i = 0; i < positions->size(); i++){
temp->insert((*((*common)[i]))[(*positions)[i]]);
}

if(temp->size() == 1) {
theIntersection->push_back(*(temp->begin()));
for(int i = 0; i < positions->size(); i++){
(*positions)[i]++;
}
}else{
auto iter = temp->end();
iter--;
int maxElement = *iter;
for(int i = 0; i < positions->size(); i++){
if((*((*common)[i]))[(*positions)[i]] < maxElement) (*positions)[i]++;
}
}

delete temp;

for(int i = 0; i < positions->size(); i++){
if((*positions)[i] < (*common)[i]->size()) continue;
else{
flag = false;
break;
}
}
}

delete positions;
return theIntersection;
}


// integer to string
string itos_c(int i){
    stringstream s;
    s << i;
    return s.str();
}


// get wall time
double GetWallTime(){
    struct timeval time;
    if(gettimeofday(&time, NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec*.000001;
}


// get CPU time
double GetCPUTime(){
    return (double)clock() / CLOCKS_PER_SEC;
}

// sort a vector of vector in descending order of vector size
bool cmp(const vector<int> &a,const vector<int> &b)
{
    return a.size() > b.size();
}

vector<bool> boolify(vector<int>& S, int n){
    vector<bool> Sbool(n, false);
    for(int i = 0; i < S.size(); i++) Sbool[S[i]] = true;
    return Sbool;
}

vector<int> FindCommon(vector<graph*> graphCollection, int i, int u){
    vector<int> vertex_list;
    vertex_list.clear();
    for(int j = 0; j < graphCollection[i]->vertices[u]->neighbors.size(); j++){
        int flag = 0;
        int cur_nb = graphCollection[i]->vertices[u]->neighbors[j];
        for(int l = 0; l < graphCollection.size(); l ++){
            if(l != i){
                if(find(graphCollection[l]->vertices[u]->kNeighbors.begin(), graphCollection[l]->vertices[u]->kNeighbors.end(), cur_nb) == graphCollection[l]->vertices[u]->kNeighbors.end()){
                    flag = 1;
                    break;
                }
            }
        }
        if(flag == 0){
            vertex_list.push_back(cur_nb);
        }
    }
    return vertex_list;
}

