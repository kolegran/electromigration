#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// Structure dimensions
#define  NX         50
#define  NY         50
#define  NZ         1

#define  Z          12                    // The number of coordination sphere

#define  E_V        1.60217656535e-19     // Electronvolt constant
#define  K_B        1.38064852e-23        // Boltzmann constant
#define  NYU_0      1.0e+13               // Debye frequency
#define  E_SADDLE   0                     // The saddle point energy
#define  A_0        0.1                   // The initial artificial noise
#define  DT         1.0e-14               // Time step
#define  C_INIT     0.5                   // Initial distribution

#define  Z_A       -10                    // The effective charge (from -10 to -30)
#define  E_X        1.0e+2                // Projection of field strength on the axis X

#define  E_CHARGE   1.60217656535e-19     // The electron charge
#define  A_SPACING  2.5e-10               // The lattice (atomic) spacing (A----B)

#define  T          550                   // The temperature

#define  V_AA      -0.32e-20              // The interaction energy

#define  STEPS      200                   // Steps for writing a file

#endif
