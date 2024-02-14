//
// Created by Hao Pan on 9/20/21.
//

#ifndef PGRAPHKCLUB_FUNCTIONS_H
#define PGRAPHKCLUB_FUNCTIONS_H


#include "classes.h"
#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <ostream>
#include <time.h>
#include <sys/time.h>
#include <string>
#include <iomanip>
#include <set>
using namespace std;

vector<string>* ReadTasks();
void ReadInput(string);
void ReadInput(int,int);
void ReadInputKClubSig(string);
vector<int> GetHeuristicSol();
graph* GetIntersectionGraph();
graph* GetIntersectionGraphOfPowerGraphs();
void Solve();
void KClubSig();
vector<int> GetPersistentKClubSig(int, int*, vector<int>*);
void CorePeel(graph*, int);
void CommunityPeel(graph*, int);
vector<vector<int>> GetComponents(graph*, int);
void DisconnectGraphsByComponents(vector<graph*>, graph*);
void DeleteEdgesBasedOnComponents(vector<graph*>, vector<int>&, vector<int>&, int&);
vector<bool> boolify(vector<int>&, int); // from parsimonious paper
vector<int> Drop(graph*, vector<bool>&); // from parsimonious paper
int KNeighborhoodSize(graph*, vector<bool>&, int); // from parsimonious paper
vector<int> FindCommon(vector<graph*>, int, int);

//vector<int> GetCore(int);
//vector<graph*> PreprocessCore();
//vector<vector<int>> GetComponents(graph*);

void CleanUp();
void CleanUpKClubSig();
//void Solve();
vector<int>* MINIMALIZE(graph*, int, int, int, vector<int>*);
vector<int>* MINIMALIZE(graph*, int, int, int, vector<int>*, vector<int>*);
bool ValidateSolution(vector<int>*);
void CalculateRatio();
bool IsAdjacentOverCollection(int, int);
bool IsKAdjacentOverCollection(int, int);
vector<int>* TakeIntersection(vector<vector<int>*>*);
string itos_c(int);
double GetWallTime();
double GetCPUTime();

bool cmp(const vector<int> &,const vector<int> &);

#endif //PGRAPHKCLUB_FUNCTIONS_H
