#include "Additional_functions.hpp"

using namespace std;

Grayatom  L[2*NX][2*NY][2*NZ];      // Lattice (argumented global variable)

// Kinetic Mean-Field method
void masterEquation(void) {
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                L[x  ][y  ][z  ].C_V_new = L[x  ][y  ][z  ].C_V + rightPartOfEquation(x,   y,   z)   * DT;
                L[x+1][y+1][z  ].C_V_new = L[x+1][y+1][z  ].C_V + rightPartOfEquation(x+1, y+1, z)   * DT;
                L[x+1][y  ][z+1].C_V_new = L[x+1][y  ][z+1].C_V + rightPartOfEquation(x+1, y,   z+1) * DT;
                L[x  ][y+1][z+1].C_V_new = L[x  ][y+1][z+1].C_V + rightPartOfEquation(x,   y+1, z+1) * DT;

                L[x  ][y  ][z  ].C_A = 1 - L[x  ][y  ][z  ].C_V_new;
                L[x+1][y+1][z  ].C_A = 1 - L[x+1][y+1][z  ].C_V_new;
                L[x+1][y  ][z+1].C_A = 1 - L[x+1][y  ][z+1].C_V_new;
                L[x  ][y+1][z+1].C_A = 1 - L[x  ][y+1][z+1].C_V_new;
            }
        }
    }

    // Lattice overwrite
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                L[x  ][y  ][z  ].C_V = L[x  ][y  ][z  ].C_V_new;
                L[x+1][y+1][z  ].C_V = L[x+1][y+1][z  ].C_V_new;
                L[x+1][y  ][z+1].C_V = L[x+1][y  ][z+1].C_V_new;
                L[x  ][y+1][z+1].C_V = L[x  ][y+1][z+1].C_V_new;
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
                    fOut << x   << "\t" << y   << "\t" << z   << "\t" << L[x  ][y  ][z  ].C_V << endl;
                    fOut << x+1 << "\t" << y+1 << "\t" << z   << "\t" << L[x+1][y+1][z  ].C_V << endl;
                    fOut << x   << "\t" << y+1 << "\t" << z+1 << "\t" << L[x  ][y+1][z+1].C_V << endl;
                    fOut << x+1 << "\t" << y   << "\t" << z+1 << "\t" << L[x+1][y  ][z+1].C_V << endl;
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
                L[x  ][y  ][z  ].C_V = C_INIT + A_0 * randomNumber();
                L[x+1][y+1][z  ].C_V = C_INIT + A_0 * randomNumber();
                L[x+1][y  ][z+1].C_V = C_INIT + A_0 * randomNumber();
                L[x  ][y+1][z+1].C_V = C_INIT + A_0 * randomNumber();
                
                L[x  ][y  ][z  ].C_A = 1 - L[x  ][y  ][z  ].C_V;
                L[x+1][y+1][z  ].C_A = 1 - L[x+1][y+1][z  ].C_V;
                L[x+1][y  ][z+1].C_A = 1 - L[x+1][y  ][z+1].C_V;
                L[x  ][y+1][z+1].C_A = 1 - L[x  ][y+1][z+1].C_V;
            }
        }   
    }
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                neighboursOfAtomPBC(x,   y,   z  );
                neighboursOfAtomPBC(x+1, y+1, z  );
                neighboursOfAtomPBC(x+1, y,   z+1);
                neighboursOfAtomPBC(x,   y+1, z+1);
            }
        }   
    }
}


// Calculation initial sum of concentrations 
void calcInitialSum(double *initSumC_A, double *initSumC_V) {
    *initSumC_A = 0.0;
    *initSumC_V = 0.0;

    for (int x = 0; x < 2*NX; x++) {
        for (int y = 0; y < 2*NY; y++) {
            for (int z = 0; z < 2*NZ; z++) {
                if ((x+y+z)%2 == 0) {
                    *initSumC_A += L[x][y][z].C_A;
                    *initSumC_V += L[x][y][z].C_V;
                }
            }
        }
    }
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

        sumOfGammaAV += (L[ni][nj][nk].C_A) * gammaAV(ai, aj, ak, ni, nj, nk);
        sumOfGammaVA +=  L[ni][nj][nk].C_V  * gammaVA(ai, aj, ak, ni, nj, nk);
    }

    return  (-L[ai][aj][ak].C_V * sumOfGammaAV + L[ai][aj][ak].C_A * sumOfGammaVA);
}


// The probability of an exchange of atoms (i, j) per unit time
double gammaAV(int ai, int aj, int ak, int ni, int nj, int nk) {
    return NYU_0 * exp(-(E_SADDLE - energyAVBefore(ai, aj, ak, ni, nj, nk) + addElectricField(ai, aj, ak, ni)) / (T*(K_B/E_V)));
}


// The probability of an exchange of atoms (j, i) per unit time
double gammaVA(int ai, int aj, int ak, int ni, int nj, int nk) {
    return NYU_0 * exp(-(E_SADDLE - energyVABefore(ai, aj, ak, ni, nj, nk) + addElectricField(ai, aj, ak, ni)) / (T*(K_B/E_V)));
}


// Sum of binding energies
double energyAVBefore(int ai, int aj, int ak, int ni, int nj, int nk) {
    return bindingEnergy(ni, nj, nk);    // E_A(In) + E_V(I), where E_V(I) = 0
}


// Sum of binding energies
double energyVABefore(int ai, int aj, int ak, int ni, int nj, int nk) {
    return bindingEnergy(ai, aj, ak);    // E_A(I) + E_V(In), where E_V(In) = 0
}


// Binding energy
double bindingEnergy(int xi, int yj, int zk) {
    int ni, nj, nk;
    double neighborsSum = 0.0;

    for (int k = 0; k < L[xi][yj][zk].numnb; k++) {
        ni            =     L[xi][yj][zk].nb[k].x;
        nj            =     L[xi][yj][zk].nb[k].y;
        nk            =     L[xi][yj][zk].nb[k].z;

        neighborsSum += 1 - L[ni][nj][nk].C_V;
    }
    return V_AA/E_V * neighborsSum;
}


// Add electric field
double addElectricField(int ax, int ay, int az, int nx) {
    return (-Z_A * E_CHARGE * E_X*(ax - nx) * A_SPACING/2.0) / E_V;
}


void compensate(void) {
    int repeat = 0;

    // Backcompenzating sites and surroundings over 1 or under 0
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                backdistribute(x,   y,   z  );
                backdistribute(x+1, y+1, z  );
                backdistribute(x+1, y,   z+1);
                backdistribute(x,   y+1, z+1);
            }
        }
    }

    // Check range [0; 1]
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                if (L[x  ][y  ][z  ].C_V < 0 || L[x  ][y  ][z  ].C_V > 1) repeat = 1;
                if (L[x+1][y+1][z  ].C_V < 0 || L[x+1][y+1][z  ].C_V > 1) repeat = 1;
                if (L[x+1][y  ][z+1].C_V < 0 || L[x+1][y  ][z+1].C_V > 1) repeat = 1;
                if (L[x  ][y+1][z+1].C_V < 0 || L[x  ][y+1][z+1].C_V > 1) repeat = 1;
            }
        }
    }

    // Backcompenzating sites and surroundings over 1 or under 0
    if (repeat == 1) {
        for (int x = 0; x < 2*NX; x += 2) {
            for (int y = 0; y < 2*NY; y += 2) {
                for (int z = 0; z < 2*NZ; z += 2) {
                    backdistribute(x,   y,   z  );
                    backdistribute(x+1, y+1, z  );
                    backdistribute(x+1, y,   z+1);
                    backdistribute(x,   y+1, z+1);
                }
            }
        }
        repeat = 0;
    }

    // Check range [0; 1]
    for (int x = 0; x < 2*NX; x += 2) {
        for (int y = 0; y < 2*NY; y += 2) {
            for (int z = 0; z < 2*NZ; z += 2) {
                if (L[x  ][y  ][z  ].C_V < 0 || L[x  ][y  ][z  ].C_V > 1) repeat = 1;
                if (L[x+1][y+1][z  ].C_V < 0 || L[x+1][y+1][z  ].C_V > 1) repeat = 1;
                if (L[x+1][y  ][z+1].C_V < 0 || L[x+1][y  ][z+1].C_V > 1) repeat = 1;
                if (L[x  ][y+1][z+1].C_V < 0 || L[x  ][y+1][z+1].C_V > 1) repeat = 1;
                if (repeat == 1) {
                    cout << "WARNING: out of range after double backdistributing" << endl;
                    cout << "C = " << L[x][y+1][z+1].C_V << endl;
                    repeat = 0;
                }
            }
        }
    }
}


// Redistribution of concentration near neighbors
void backdistribute(int x, int y, int z) {
    double sum, diff, d_comp;
    int ni, nj, nk;

    if (L[x][y][z].C_V < 0.0) {
        sum = 0;
        for (int k = 0; k < Z ; k++) {
            ni = L[x][y][z].nb[k].x;
            nj = L[x][y][z].nb[k].y;
            nk = L[x][y][z].nb[k].z;
            sum += L[ni][nj][nk].C_V;
        }

        diff = L[x][y][z].C_V;
        L[x][y][z].C_V = 0.0;
        L[x][y][z].C_A = 1.0;

        for (int k = 0; k < Z ; k++) {
            ni = L[x][y][z].nb[k].x;
            nj = L[x][y][z].nb[k].y;
            nk = L[x][y][z].nb[k].z;

            d_comp = diff * L[ni][nj][nk].C_V/sum;
            L[ni][nj][nk].C_V += d_comp;
            L[ni][nj][nk].C_A = 1 - L[ni][nj][nk].C_V;
        }

    } else if (L[x][y][z].C_V > 1.0) {
        sum = 0;
        for (int k = 0; k < Z ; k++) {
            ni = L[x][y][z].nb[k].x;
            nj = L[x][y][z].nb[k].y;
            nk = L[x][y][z].nb[k].z;
            sum += L[ni][nj][nk].C_V;
        }

        diff = L[x][y][z].C_V - 1.0;
        L[x][y][z].C_V = 1.0;
        L[x][y][z].C_A = 0.0;

        for (int k = 0; k < Z ; k++) {
            ni = L[x][y][z].nb[k].x;
            nj = L[x][y][z].nb[k].y;
            nk = L[x][y][z].nb[k].z;

            d_comp = diff * (1 - L[ni][nj][nk].C_V)/(Z - sum);
            L[ni][nj][nk].C_V += d_comp;
            L[ni][nj][nk].C_A = 1 - L[ni][nj][nk].C_V;
        }
    }
}


// Сheck for substance retention
void checkConservationLaws(double initialSumC_A, double initialSumC_V) {
    double epsilon1 = 1e-1;
    double epsilon2 = 1e-1;
    double currentSumC_A = 0.0;     // Current sum of concentration C_A
    double currentSumC_V = 0.0;     // Current sum of concentration C_V
    
    
    // Calculate current sum of concentrations (C_A, C_B)
    for (int x = 0; x < 2*NX; x++) {
        for (int y = 0; y < 2*NY; y++) {
            for (int z = 0; z < 2*NZ; z++) {
                if ( (x+y+z)%2 == 0 ) {
                    currentSumC_A += L[x][y][z].C_A;
                    currentSumC_V += L[x][y][z].C_V;
                }
            }
        }
    }

    // Sum of concentrations C_A must be const
    if (fabs(initialSumC_A - currentSumC_A) > epsilon1) {
        cout << "Error: Sum of concentrations C_A are not constant!" << endl;
        cout << "initialSumC_A = " << initialSumC_A << endl;
        cout << "currentSumC_A = " << currentSumC_A << endl;
        cin.get();
        exit(0);
    }

    // Sum of concentrations C_V must be const
    if (fabs(initialSumC_V - currentSumC_V) > epsilon1) {
        cout << "Error: Sum of concentrations C_V are not constant!" << endl;
        cout << "initialSumC_V = " << initialSumC_V << endl;
        cout << "currentSumC_V = " << currentSumC_V << endl;
        cin.get();
        exit(0);
    }
    
    // 0 > concentration > 1
    for (int x = 0; x < 2*NX; x++) {
        for (int y = 0; y < 2*NY; y++) {
            for (int z = 0; z < 2*NZ; z++) {
                if ((x+y+z)%2 == 0) {
                    // C_A < 0
                    if (L[x][y][z].C_A < 0.0 - epsilon2) {
                        cout << "Error: " << x << ", " << y << ", " << z << ", " << L[x][y][z].C_A << endl;
                        cout << "Concentration C_A < 0" << endl;
                        cin.get();
                        exit(0);
                    }
                    // C_A > 1
                    if (L[x][y][z].C_A > 1.0 + epsilon2) {
                        cout << "Error: " << x << ", " << y << ", " << z << ", " << L[x][y][z].C_A << endl;
                        cout << "Concentration C_A > 1" << endl;
                        cin.get();
                        exit(0);
                    }
                    // C_V < 0
                    if (L[x][y][z].C_V < 0.0 - epsilon2) {
                        cout << "Error: " << x << ", " << y << ", " << z << ", " << L[x][y][z].C_V << endl;
                        cout << "Concentration C_V < 0" << endl;
                        cin.get();
                        exit(0);
                    }
                    // C_V > 1
                    if (L[x][y][z].C_V > 1.0 + epsilon2) {
                        cout << "Error: " << x << ", " << y << ", " << z << ", " << L[x][y][z].C_V << endl;
                        cout << "Concentration C_V > 1" << endl;
                        cin.get();
                        exit(0);
                    }
                }
            }
        }
    }
}
