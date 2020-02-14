// this will be the header file
// creates initial conditions of random ('r'), Gaussian ('n'), or a mixture of the Gaussians ('m')

#ifndef IC_GENERATOR
#define IC_GENERATOR

#include <iostream>
#include <lib/AliTools/Psi_Zeta_Init_Class.h>

//#include <fstream>
//#include <list>

#include <lib/arrays/ArrayV.h>
#include <lib/amr/OcTree.h>
#include <lib/amr/OcTreeCellBased.h>
//#include <lib/amr/OcTreeLevelSet.h>
//#include <lib/amr/OcTreeAdvection.h>
//#include <lib/tools/OcTreeToolbox.h>
//#include <lib/amr/OcTreeCellBasedLevelSet.h>
//#include <lib/solvers/OcTreeSolverCellCenteredPoisson.h>

//#include <sys/stat.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>

//#include <stdio.h>
//#include <sys/uio.h>

//#include </Users/aliheydari/stdc++.h>
//#include <iostream>
//#include <sys/stat.h>
//#include <sys/types.h>
//#include<unistd.h>

//#include <mains_prions/prion_cell_splitting.cpp>



using namespace std;
using namespace CASL;

// Structs for weight ditribution
struct IC_Generator {
    ArrayV<double> zeta, psi;

//    IC_Generator generator(char choice ,ArrayV<double> current_A,ArrayV<double> current_B,OcTreeCellBased my_octree, Initial_Solution psi_initial, zeta_Initial_Solution zeta_initial);

};


IC_Generator gen(char choice ,ArrayV<double> current_A,ArrayV<double> current_B,OcTreeCellBased my_octree, Initial_Solution psi_initial, zeta_Initial_Solution zeta_initial);

#endif // IC_GENERATOR
