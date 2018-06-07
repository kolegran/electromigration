#ifndef ADDITIONAL_FUNCTIONS_HPP
#define ADDITIONAL_FUNCTIONS_HPP

#include <iostream>                                     // cout(), cin()
#include <iomanip>                                      // setfill(), setw()
#include <cstdlib>                                      // rand(), srand()
#include <ctime>                                        // time()
#include <cmath>                                        // exp(), fabs()
#include <fstream>                                      // open(), close()
#include <string>
#include <sstream>                                      // str()

#include "Constants.hpp"

// Coordinates of a neighboring atom
typedef struct Nbatom {
    int x;
    int y;
    int z;
} Nbatom;


// Gray atom
typedef struct GRAYATOM {
    double C_old, C_new;    // Old and new concentrations
    Nbatom nb[Z];           // Coordinates of neighboring atoms
    int numnb;              // The number of neighbors
} Grayatom;

                                                                        // ---function prototypes
double randomNumber(void);                                              // Generation of a random number [0; 1]
void   neighboursOfAtomPBC(int x, int y, int z);                        // Memorizing the neighbors of each atom with periodic boundary conditions
void   latticeInitialization(void);                                     // Lattice initialization

void   masterEquation(void);
double rightPartOfEquation(int ai, int aj, int ak);                     // The right-hand side of equation dCi/dt = ...

double gammaAV(int ai, int aj, int ak, int ni, int nj, int nk);         // The probability of an exchange of atoms (i, j) per unit time
double gammaVA(int ai, int aj, int ak, int ni, int nj, int nk);         // The probability of an exchange of atoms (i, j) per unit time

double energyAVBefore(int ai, int aj, int ak, int ni, int nj, int nk);  // Sum of binding energies
double energyVABefore(int ai, int aj, int ak, int ni, int nj, int nk);  // Sum of binding energies
double bindingEnergy(int xi, int yj, int zk);                           // Binding energy A-B (B-A)

double addElectricField(int ax, int nx);                                // Add electric field

double calcInitialSum(void);                         
void   checkConservationLaws(double initialSum);                        // For check conservation laws

void   fileXYZ(void);                                                   // write data to the file

#endif

