#include "Additional_functions.hpp"

using namespace std;

int main(void) {
    double endTime = 10.0;    // End time
    double initSum = 0.0;  
    int    stepCtn = 0;       // Counter for steps

    srand(time(nullptr));

    latticeInitialization();

    initSum = calcInitialSum();

    fileXYZ();
    
    // The main loop
    for (double time = 0; time < endTime; time += DT) {
        masterEquation();
        
        checkConservationLaws(initSum);
        
        // Writing to a file
        if (!(stepCtn%STEPS)) {
            fileXYZ();
        }
        stepCtn++;
    }
    cin.get();
    return 0;
}
