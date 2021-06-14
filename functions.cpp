//
// Created by Hao Pan on 4/5/21.
//

#include "functions.h"
#include "classes.h"

//Color::Modifier green(Color::FG_GREEN);
//Color::Modifier def(Color::FG_DEFAULT);


vector<string>* ReadTasks(){
    
    vector<string>* tasks = new vector<string>();
    string temp;
    ifstream fin("tasks.txt");
    while(getline(fin, temp))
        tasks->push_back(temp);
    return tasks;
}


vector<graph*> graphCollection;
int pParam, kParam, firstGraph, lastGraph;
string instanceName, methodParam;
void ReadInput(string task){

    stringstream ss(task);
    int u, v;
    string temp;

    ss >> pParam >> firstGraph >> lastGraph >> instanceName >> kParam >> methodParam;
    //cout << green << "\nREAD IN INPUT..." << def << endl;
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


void Solve(){

    //cout << green << "\nSOLVE FOR A MAX " << pParam << "-GRAPH " << kParam << "-CLUB... " << def << endl;
    cout << "\nSOLVE FOR A MAX " << pParam << "-GRAPH " << kParam << "-CLUB... " << endl;

    try{

        double timeBegin = GetWallTime();

        for(auto graphPtr : graphCollection)
            graphPtr->FindKNeighbors(kParam);

        GRBEnv* ENV = 0;
        GRBVar* X = 0;

        ENV = new GRBEnv();
        GRBModel MODEL = GRBModel(*ENV);
        int nParam = graphCollection[0]->n;

        X = MODEL.addVars(nParam, GRB_BINARY);
        MODEL.update();

        for(int i = 0; i < nParam; i++)
            X[i].set(GRB_DoubleAttr_Obj, 1.0);
        MODEL.update();

        for(int i = 0; i < nParam; i++){
            for(int j = i + 1; j < nParam; j++){
                int goFlag = 0;
                for(auto graphPtr : graphCollection){
                    if(graphPtr->IsKAdjacent(i, j) == false){
                        goFlag = 1;
                        break;
                    }
                }
                if(goFlag == 1)
                    MODEL.addConstr(X[i] + X[j] <= 1);
            }
        }
        MODEL.update();

        MODEL.set(GRB_IntParam_MIPFocus, 0);
        MODEL.set(GRB_IntParam_Threads, 0);
        MODEL.set(GRB_IntParam_Method, -1);
        MODEL.set(GRB_IntParam_NodeMethod, 1);
        MODEL.set(GRB_IntParam_Cuts, 0);
        MODEL.set(GRB_DoubleParam_TimeLimit, 7200);
        MODEL.set(GRB_DoubleParam_MIPGap, 0);
        MODEL.set(GRB_StringParam_LogFile, "");
        MODEL.set(GRB_IntParam_OutputFlag, 1);
        MODEL.set(GRB_IntAttr_ModelSense, -1);
        MODEL.set(GRB_IntParam_LazyConstraints,1);
        int numStrengthenedLazyCut = 0, numLazyCut = 0, numCallback = 0;
        double callbackTime = 0, aveNumTermsNeg = 0, aveNumTermsNegStrengthened = 0;
        pGraph_callback myCallback = pGraph_callback(methodParam, X, graphCollection, kParam, numStrengthenedLazyCut, numLazyCut, numCallback, callbackTime, aveNumTermsNeg, aveNumTermsNegStrengthened);
        MODEL.setCallback(&myCallback);

        MODEL.optimize();

        double timeEnd = GetWallTime();

        //validate solution
        if(MODEL.get(GRB_IntAttr_SolCount)){
            vector<int>* solutionPtr = new vector<int>();
            for(int i = 0; i < nParam; i++){
                if(X[i].get(GRB_DoubleAttr_X) > 0.5)
                    solutionPtr->push_back(i);
            }

            if(ValidateSolution(solutionPtr) == false){
                delete solutionPtr;
                return;
            }
            delete solutionPtr;
        }else{
            cout << "NO SOLUTION FOUND SO FAR" << endl;
            return;
        }

        //cout and write out results
        double objValue = MODEL.get(GRB_DoubleAttr_ObjVal);
        double objBound = MODEL.get(GRB_DoubleAttr_ObjBound);
        int numBBNode = MODEL.get(GRB_DoubleAttr_NodeCount);
        double duration = timeEnd - timeBegin;

        delete ENV;

        //cout << green << "\nSUMMARY" << def << endl;
        cout << "\nSUMMARY" << endl;
        if(methodParam == "BASE"){

            cout << "PARAMETER P : " << pParam << endl;
            cout << "PARAMETER K : " << kParam << endl;
            cout << "INSTANCE : " << instanceName << " " << firstGraph << "-" << lastGraph << endl;
            cout << "METHOD : " << methodParam << endl;
            cout << "OBJ VALUE : " << objValue << endl;
            cout << "OBJ BOUND : " << objBound << endl;
            cout << "DURATION : " << fixed << setprecision(2) << duration << endl;
            cout << "NUM OF CALLBACKS : " << numCallback << endl;
            cout << "CALLBACK TIME : " << fixed << setprecision(2) << callbackTime << endl;
            cout << "NUM OF LAZY CONSTRAINTS ADDED : " << numLazyCut << endl;
            cout << "AVE NUM OF NEGATIVE TERMS : " << aveNumTermsNeg << endl;
            cout << "NUM OF BB NODES : " << numBBNode << endl;

            ofstream fout;
            fout.open("BASE.csv", ios::in);
            if(fout.fail()){
                fout.open("BASE.csv", ios::out);
                fout << "p, k, instance, obj, obj bound, duration, num callbacks, callback time, num lazy cuts, ave num of negative terms, num BB nodes\n";
                fout.close();
            }else{
                fout.close();
            }

            fout.open("BASE.csv", ios::app);
            fout << pParam << ",";
            fout << kParam << ",";
            fout << instanceName << " " << firstGraph << "-" << lastGraph << ",";
            fout << objValue << ",";
            fout << objBound << ",";
            fout << fixed << setprecision(2) << duration << ",";
            fout << numCallback << ",";
            fout << fixed << setprecision(2) << callbackTime << ",";
            fout << numLazyCut << ",";
            fout << aveNumTermsNeg << ",";
            fout << numBBNode << endl;
            fout.close();
        }else{

            cout << "PARAMETER P : " << pParam << endl;
            cout << "PARAMETER K : " << kParam << endl;
            cout << "INSTANCE : " << instanceName << " " << firstGraph << "-" << lastGraph << endl;
            cout << "METHOD : " << methodParam << endl;
            cout << "OBJ VALUE : " << objValue << endl;
            cout << "OBJ BOUND : " << objBound << endl;
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
                fout << "p, k, instance, obj, obj bound, duration, num callbacks, callback time, num lazy cuts, num strengthened lazy cuts, ave num of negative terms, ave num of negative terms strengthened, num BB nodes\n";
                fout.close();
            }else{
                fout.close();
            }

            fout.open(methodParam + ".csv", ios::app);
            fout << pParam << ",";
            fout << kParam << ",";
            fout << instanceName << " " << firstGraph << "-" << lastGraph << ",";
            fout << objValue << ",";
            fout << objBound << ",";
            fout << fixed << setprecision(2) << duration << ",";
            fout << numCallback << ",";
            fout << fixed << setprecision(2) << callbackTime << ",";
            fout << numLazyCut << ",";
            fout <<  numStrengthenedLazyCut << ",";
            fout << aveNumTermsNeg << ",";
            fout << aveNumTermsNegStrengthened << ",";
            fout << numBBNode << endl;
            fout.close();
        }
        //cout << green << "\nFINISH WRITING OUT RESULTS" << def << endl;
        cout << "\nFINISH WRITING OUT RESULTS" << endl;
    }catch (GRBException e) {
        cout << "ERROR CODE = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "EXCEPTION DURING OPTIMIZATION" << endl;
    }
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
