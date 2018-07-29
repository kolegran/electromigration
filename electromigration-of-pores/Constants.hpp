#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// Structure dimensions
#define  NX         20
#define  NY         20
#define  NZ         15

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

#define  U          1.0e+8*A_SPACING*NX   // U ~ E*d ~ 1e+8 * NX * a ~ 2.5e-2 * NX
#define  C_SIGMA    1.0e+7                // Electrical conductivity (c - constant)
#define  PHI_LEFT  -U/2.0
#define  PHI_RIGHT  U/2.0
#define  EPS        1.0e-3

#define  T          550                   // The temperature

#define  V_AA      -0.32e-20              // The interaction energy

#define  STEPS      350                   // Steps for writing a file

#endif
