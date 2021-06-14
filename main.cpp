

#include "functions.h"

int main() {

    vector<string>* tasks = ReadTasks();
    for(auto task : (*tasks)){
        ReadInput(task);
        Solve();
        CleanUp();
    }
    delete tasks;


//    vector<string>* tasks = ReadTasks();
//    for(auto task : (*tasks)){
//        ReadInput(task);
//        CalculateRatio();
//        CleanUp();
//    }
//    delete tasks;

    return 0;
}
