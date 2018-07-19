#include "Additional_functions.hpp"

using namespace std;

int main(void) {
    double endTime = 10.0;       // End time 
    double initialSumC_A = 0.0;
    double initialSumC_V = 0.0;
    int    stepCtn = 0;          // Counter for steps

    srand(time(nullptr));

    latticeInitialization();

    //calcInitialSum(&initialSumC_A, &initialSumC_V);

    fileXYZ();
    
    // The main loop
    for (double time = 0; time < endTime; time += DT) {
        masterEquation();
        
        compensate();

        //checkConservationLaws(initialSumC_A, initialSumC_V);
        
        // Writing to a file
        if (!(stepCtn%STEPS)) {
            fileXYZ();
        }
        stepCtn++;
    }
    cin.get();
    return 0;
}
