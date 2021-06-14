//
// Created by Hao Pan on 4/5/21.
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
void CleanUp();
void Solve();
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


//namespace Color{
//    enum Code {
//        FG_RED      = 31,
//        FG_GREEN    = 32,
//        FG_BLUE     = 34,
//        FG_DEFAULT  = 39,
//        BG_RED      = 41,
//        BG_GREEN    = 42,
//        BG_BLUE     = 44,
//        BG_DEFAULT  = 49
//    };
//    class Modifier {
//        Code code;
//    public:
//        Modifier(Code pCode) : code(pCode) {}
//        friend std::ostream&
//        operator<<(std::ostream& os, const Modifier& mod) {
//            return os << "\033[" << mod.code << "m";
//        }
//    };
//}

#endif //PGRAPHKCLUB_FUNCTIONS_H
