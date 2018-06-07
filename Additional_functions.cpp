#include "Additional_functions.hpp"

using namespace std;

Grayatom  L[2*NX][2*NY][2*NZ];      // Lattice (argumented global variable)

// Kinetic Mean-Field method
void masterEquation(void) {
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                L[x  ][y  ][z  ].C_new = L[x  ][y  ][z  ].C_old + rightPartOfEquation(x,   y,   z)   * DT;
                L[x+1][y+1][z  ].C_new = L[x+1][y+1][z  ].C_old + rightPartOfEquation(x+1, y+1, z)   * DT;
                L[x+1][y  ][z+1].C_new = L[x+1][y  ][z+1].C_old + rightPartOfEquation(x+1, y,   z+1) * DT;
                L[x  ][y+1][z+1].C_new = L[x  ][y+1][z+1].C_old + rightPartOfEquation(x,   y+1, z+1) * DT;
            }
        }
    }

    // Lattice overwrite
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                L[x  ][y  ][z  ].C_old = L[x  ][y  ][z  ].C_new;
                L[x+1][y+1][z  ].C_old = L[x+1][y+1][z  ].C_new;
                L[x+1][y  ][z+1].C_old = L[x+1][y  ][z+1].C_new;
                L[x  ][y+1][z+1].C_old = L[x  ][y+1][z+1].C_new;
            }
        }
    }
}


// Write data to the file
void fileXYZ(void) {
    static int fnameCtn = 0;    // Counter for names of files

    fnameCtn++;
    // Create name of the file
    ostringstream s(ostringstream::ate);
    s << setfill('0') << setw(9) << fnameCtn << ".txt";

    // Create the file
    ofstream fOut;

    // Open the file
    fOut.open(s.str());
 
    if (fOut) {
        // Write in file the number of atoms
        fOut << 4*NX*NY*NZ << endl << endl;

        for (int x = 0; x < 2*NX; x += 2) {
            for (int y = 0; y < 2*NY; y += 2) {
                for (int z = 0; z < 2*NZ; z += 2) {
                    // Write in file coordinates of atoms
                    fOut << x   << "\t" << y   << "\t" << z   << "\t" << L[x  ][y  ][z  ].C_old << endl;
                    fOut << x+1 << "\t" << y+1 << "\t" << z   << "\t" << L[x+1][y+1][z  ].C_old << endl;
                    fOut << x   << "\t" << y+1 << "\t" << z+1 << "\t" << L[x  ][y+1][z+1].C_old << endl;
                    fOut << x+1 << "\t" << y   << "\t" << z+1 << "\t" << L[x+1][y  ][z+1].C_old << endl;
                }
            }
        }
    } else {
        cout << "Warning: couldn\'t open the file" << endl;
    }
    // Close the file
    fOut.close();
}


// Generation of a random number [0; 1]
double randomNumber(void) {
    return 2*(double)rand() / (double)RAND_MAX - 1;
}


// Lattice initialization
// Check periodic boundary conditions
void latticeInitialization(void) {
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                L[x  ][y  ][z  ].C_old = C_INIT + A_0 * randomNumber();
                neighboursOfAtomPBC(x,   y,   z);

                L[x+1][y+1][z  ].C_old = C_INIT + A_0 * randomNumber();
                neighboursOfAtomPBC(x+1, y+1, z);

                L[x+1][y  ][z+1].C_old = C_INIT + A_0 * randomNumber();
                neighboursOfAtomPBC(x+1, y,   z+1);

                L[x  ][y+1][z+1].C_old = C_INIT + A_0 * randomNumber();
                neighboursOfAtomPBC(x,   y+1, z+1);
            }
        }   
    }
}


// Calculation initial sum of concentrations 
double calcInitialSum(void) {
    double concentrationInitSum = 0.0;
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                concentrationInitSum += L[x  ][y  ][z  ].C_old;
                concentrationInitSum += L[x+1][y+1][z  ].C_old;
                concentrationInitSum += L[x+1][y  ][z+1].C_old;
                concentrationInitSum += L[x  ][y+1][z+1].C_old;
            }
        }
    }
    return concentrationInitSum;
}


// Memorizing the neighbors of each atom
void neighboursOfAtomPBC(int x, int y, int z) {
    int xm1, xp1;                // m-minus, p-plus
    int ym1, yp1;                // ym1 = y-1; yp1 = y+1
    int zm1, zp1;                // zm1 = z-1; zp1 = z+1
    int i = 0;                   // The counter of the neighbors

    if (x == 2*NX-1)  xp1 = 0;  else  xp1 = x + 1;
    if (y == 2*NY-1)  yp1 = 0;  else  yp1 = y + 1;
    if (z == 2*NZ-1)  zp1 = 0;  else  zp1 = z + 1;

    if (x == 0)  xm1 = 2*NX-1;  else  xm1 = x - 1;
    if (y == 0)  ym1 = 2*NY-1;  else  ym1 = y - 1;
    if (z == 0)  zm1 = 2*NZ-1;  else  zm1 = z - 1;

    L[x][y][z].nb[i].x = xm1;  L[x][y][z].nb[i].y = ym1;  L[x][y][z].nb[i].z = z;   i++;
    L[x][y][z].nb[i].x = xm1;  L[x][y][z].nb[i].y = yp1;  L[x][y][z].nb[i].z = z;   i++;
    L[x][y][z].nb[i].x = xp1;  L[x][y][z].nb[i].y = ym1;  L[x][y][z].nb[i].z = z;   i++;
    L[x][y][z].nb[i].x = xp1;  L[x][y][z].nb[i].y = yp1;  L[x][y][z].nb[i].z = z;   i++;

    L[x][y][z].nb[i].x = x;    L[x][y][z].nb[i].y = ym1;  L[x][y][z].nb[i].z = zm1; i++;
    L[x][y][z].nb[i].x = x;    L[x][y][z].nb[i].y = ym1;  L[x][y][z].nb[i].z = zp1; i++;
    L[x][y][z].nb[i].x = x;    L[x][y][z].nb[i].y = yp1;  L[x][y][z].nb[i].z = zm1; i++;
    L[x][y][z].nb[i].x = x;    L[x][y][z].nb[i].y = yp1;  L[x][y][z].nb[i].z = zp1; i++;

    L[x][y][z].nb[i].x = xm1;  L[x][y][z].nb[i].y = y;    L[x][y][z].nb[i].z = zm1; i++;
    L[x][y][z].nb[i].x = xm1;  L[x][y][z].nb[i].y = y;    L[x][y][z].nb[i].z = zp1; i++;
    L[x][y][z].nb[i].x = xp1;  L[x][y][z].nb[i].y = y;    L[x][y][z].nb[i].z = zm1; i++;
    L[x][y][z].nb[i].x = xp1;  L[x][y][z].nb[i].y = y;    L[x][y][z].nb[i].z = zp1; i++;

    L[x][y][z].numnb = i;          // The number of the neighbors
}


// The right-hand side of equation dCi/dt = ...
double rightPartOfEquation(int ai, int aj, int ak) {
    int ni, nj, nk; 
    double sumOfGammaAV = 0.0;
    double sumOfGammaVA = 0.0;

    for (int k = 0; k < L[ai][aj][ak].numnb; k++) {
        ni = L[ai][aj][ak].nb[k].x;
        nj = L[ai][aj][ak].nb[k].y;
        nk = L[ai][aj][ak].nb[k].z;
        sumOfGammaAV += (1 - L[ni][nj][nk].C_old) * gammaAV(ai, aj, ak, ni, nj, nk);
        sumOfGammaVA +=     L[ni][nj][nk].C_old   * gammaVA(ai, aj, ak, ni, nj, nk);
    }

    return  - L[ai][aj][ak].C_old        * sumOfGammaAV
            + (1 - L[ai][aj][ak].C_old)  * sumOfGammaVA;
}


// The probability of an exchange of atoms (i, j) per unit time
double gammaAV(int ai, int aj, int ak, int ni, int nj, int nk) {
    return NYU_0 * exp(-(E_SADDLE - energyAVBefore(ai, aj, ak, ni, nj, nk) + addElectricField(ai, ni)) / THETA);
}


// The probability of an exchange of atoms (j, i) per unit time
double gammaVA(int ai, int aj, int ak, int ni, int nj, int nk) {
    return NYU_0 * exp(-(E_SADDLE - energyVABefore(ai, aj, ak, ni, nj, nk) + addElectricField(ai, ni)) / THETA);
}


// Sum of binding energies
double energyAVBefore(int ai, int aj, int ak, int ni, int nj, int nk) {
    return bindingEnergy(ni, nj, nk);
}


// Sum of binding energies
double energyVABefore(int ai, int aj, int ak, int ni, int nj, int nk) {
    return bindingEnergy(ai, aj, ak);
}


// Binding energy
double bindingEnergy(int xi, int yj, int zk) {
    int ni, nj, nk;
    double neighborsSum = 0.0;

    for (int k = 0; k < L[xi][yj][zk].numnb; k++) {
        ni            =     L[xi][yj][zk].nb[k].x;
        nj            =     L[xi][yj][zk].nb[k].y;
        nk            =     L[xi][yj][zk].nb[k].z;
        neighborsSum += 1 - L[ni][nj][nk].C_old;
    }
    return V_AA * neighborsSum/THETA;
}


// Add electric field
double addElectricField(int ax, int nx) {
    return -Z_A * E_CHARGE * E_X*(ax - nx) * A_SPACING/2.0;
}


// Сheck for substance retention
void checkConservationLaws(double initialSum) {
    double epsilon1 = 5e-1;
    double epsilon2 = 1e-1;
    double concentrationCurrentSum = 0.0;     // Current sum of concentrations

    // Calculate current sum of concentrations
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                concentrationCurrentSum += L[x  ][y  ][z  ].C_old;
                concentrationCurrentSum += L[x+1][y+1][z  ].C_old;
                concentrationCurrentSum += L[x+1][y  ][z+1].C_old;
                concentrationCurrentSum += L[x  ][y+1][z+1].C_old;
            }
        }
    }

    // Sum of concentrations must be const
    if (fabs(initialSum - concentrationCurrentSum) > epsilon1) {
        cout << "Error: Sum of concentrations are not constant!" << endl;
        cout << "concentrationInitSum = "    << initialSum << endl;
        cout << "concentrationCurrentSum = " << concentrationCurrentSum << endl;
        cin.get();
        exit(0);
    }


    // 0 < concentration < 1
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                // ... < 0
                if (L[x  ][y  ][z  ].C_old  < 0.0 - epsilon2) {
                    cout << "Error: " << x   << ", " << y   << ", " << z   << ", " << L[x  ][y  ][z  ].C_old << endl;
                    cout << "Concentration < 0" << endl;
                    cin.get();
                    exit(0);
                }
                if (L[x+1][y+1][z  ].C_old  < 0.0 - epsilon2) {
                    cout << "Error: " << x+1 << ", " << y+1 << ", " << z   << ", " << L[x+1][y+1][z  ].C_old << endl;
                    cout << "Concentration < 0" << endl;
                    cin.get();
                    exit(0);
                }
                if (L[x+1][y  ][z+1].C_old  < 0.0 - epsilon2) {
                    cout << "Error: " << x+1 << ", " << y   << ", " << z+1 << ", " << L[x+1][y  ][z+1].C_old << endl;
                    cout << "Concentration < 0" << endl;
                    cin.get();
                    exit(0);
                }
                if (L[x  ][y+1][z+1].C_old  < 0.0 - epsilon2) {
                    cout << "Error: " << x   << ", " << y+1 << ", " << z+1 << ", " << L[x  ][y+1][z+1].C_old << endl;
                    cout << "Concentration < 0" << endl;
                    cin.get();
                    exit(0);
                }
                // ... > 1
                if (L[x  ][y  ][z  ].C_old  > 1.0 + epsilon2) {
                    cout << "Error: " << x   << ", " << y   << ", " << z   << ", " << L[x  ][y  ][z  ].C_old << endl;
                    cout << "Concentration > 1" << endl;
                    cin.get();
                    exit(0);
                }
                if (L[x+1][y+1][z  ].C_old  > 1.0 + epsilon2) {
                    cout << "Error: " << x+1 << ", " << y+1 << ", " << z   << ", " << L[x+1][y+1][z  ].C_old << endl;
                    cout << "Concentration > 1" << endl;
                    cin.get();
                    exit(0);
                }
                if (L[x+1][y  ][z+1].C_old  > 1.0 + epsilon2) {
                    cout << "Error: " << x+1 << ", " << y   << ", " << z+1 << ", " << L[x+1][y  ][z+1].C_old << endl;
                    cout << "Concentration > 1" << endl;
                    cin.get();
                    exit(0);
                }
                if (L[x  ][y+1][z+1].C_old  > 1.0 + epsilon2) {
                    cout << "Error: " << x   << ", " << y+1 << ", " << z+1 << ", " << L[x  ][y+1][z+1].C_old << endl;
                    cout << "Concentration > 1" << endl;
                    cin.get();
                    exit(0);
                }
            }
        }
    }
}
