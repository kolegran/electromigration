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
    double C_A, C_A_new;
    double C_V, C_V_new;
    double phiOld, phiNew;
    Nbatom nb[Z];           // Coordinates of neighboring atoms
    int numnb;              // The number of neighbors
} Grayatom;

double randomNumber(void);                                                  // Generation of a random number (0; 1)
void   neighborsOfAtomPBC(int x, int y, int z);                             // Memorizing the neighbors of each atom with periodic boundary conditions
void   zeroConcentrations(void);                                            // Zeroing of concentrations on the faces
void   potentialsInitialization(void);                                      // Potentials initialization
void   latticeInitialization(void);                                         // Lattice initialization

void   masterEquation(void);
double rightPartOfEquation(int ai, int aj, int ak);                         // The right-hand side of equation dCi/dt = ...

double gammaAV(int ai, int aj, int ak, int ni, int nj, int nk);             // The probability of an exchange of atoms (i, j) per unit time
double gammaVA(int ai, int aj, int ak, int ni, int nj, int nk);             // The probability of an exchange of atoms (i, j) per unit time

double energyAVBefore(int ai, int aj, int ak, int ni, int nj, int nk);      // Sum of binding energies
double energyVABefore(int ai, int aj, int ak, int ni, int nj, int nk);      // Sum of binding energies
double bindingEnergy(int xi, int yj, int zk);                               // Binding energy A-B (B-A)

void   calcInitialSum(double *initSumC_A, double *initSumC_V);
void   checkConservationLaws(double initialSumC_A, double initialSum_C_V);  // For check conservation laws

void   fileXYZ(void);                                                       // Write data to the file
void   initFromFile(void);                                                  // Initialization from the file

void   compensate(void);                                                    // Recompensation of concentration in the case of out of range [0; 1]
void   backdistribute(int x, int y, int z);                                 // Backdistribute of concentration in the case of out of range [0; 1]

double addElectricField(int ax, int ay, int az, int nx, int ny, int nz);    // Add electric field
double calculationSigma(int ax, int ay, int az);
bool   checkSumCondition(void);
void   calculationPotential(void);

#endif
