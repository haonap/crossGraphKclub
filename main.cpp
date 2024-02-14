#include "functions.h"

int main() {
    vector<string>* tasks = ReadTasks();
    for(auto task : (*tasks)){

        //for cross Graph KClub
        ReadInput(task);
        Solve();
        CleanUp();

        /*
        //for KClubSig to generate Table 5-7
        ReadInputKClubSig(task);
        KClubSig();
        CleanUpKClubSig();
         */

}
delete tasks;

return 0;
}