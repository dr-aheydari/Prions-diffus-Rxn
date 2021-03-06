#ifdef CASL_OPENMP
#include <omp.h>
#define TRUE  1
#define FALSE 0
#else
#define omp_get_thread_num() 0
#endif

#include <iostream>
#include <fstream>
#include <list>
#include <tuple>

#include <lib/arrays/ArrayV.h>
#include <lib/amr/OcTree.h>
#include <lib/amr/OcTreeCellBased.h>
#include <lib/amr/OcTreeLevelSet.h>
#include <lib/amr/OcTreeAdvection.h>
#include <lib/tools/OcTreeToolbox.h>
#include <lib/amr/OcTreeCellBasedLevelSet.h>
#include <lib/solvers/OcTreeSolverCellCenteredPoisson.h>

// My custom tools
#include <lib/AliTools/isNan.h>
#include <lib/AliTools/Make_Folders.h>
#include <lib/AliTools/IC_Generator.h>
#include <lib/AliTools/Psi_Zeta_Init_Class.cpp>
#include <lib/AliTools/signum.cpp>

#ifdef CASL_PETSC
#include <lib/tools/PetscTools.h>
#endif

#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <stdio.h>
#include <sys/uio.h>

#include </Users/aliheydari/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include<unistd.h>


using namespace std;
using namespace CASL;

OcTreeLevelSet test;

OcTreeCellBased my_octree; //gives a coordinate system
OcTreeCellBased new_octree;


//array defined by Maxime
ArrayV <double> psi_n; //what we're solving for
ArrayV <double> psi_nm1; //

//ZETA VARIABLES
ArrayV <double> zeta_n;
ArrayV <double> zeta_nm1;


// Time discretization
double T = 9.5; //final time
double dt = T/500.; //time step
double tOrder = dt;
double t = 0; //kind of the same as dt
int n = 0; // just for iterations

//// Space discretization
double xL = 0.75; //space [0,L]
double yL = 0.75;
double zL = 0.75;


////////////
// run min 5, max 7 for camera ready

int min_level = 2; //for our adaptive meshing
int max_level = 6; //for our adaptive meshing
////////////

// System parameters
double D_psi = 10; // Diffusion coeffcient for "healthy" protein
double D_zeta = 0.01; // Diffusion coeffcient for aggregate

//rate
double gamma_AtoB = 0.1 ;
double initial_pop = 10*0;
double gamma_BtoA = 0.01;
double mu = 0.2 * 0;
//double D_z = 0.001;

double init_mass_A;
double init_mass_B;

double epsilon = 0.000000001;

// this threshold is used for knowing when the daughter cell will leave the mother cell
double threshold = 9;

// intialize for explicit task dependent file IO
char hole = 'T';
char DiffRXN = 'F';
char DiffOnly='F';
char compart = 'T';

char* FolderPath;
string txtPathA;
string txtPathB;
string FullPath;
string daughterPathA;
string daughterPathB;
string motherPathA;
string motherPathB;
string Return_FolderPAth;

// here we set a variable that will tell a struct which IC we have:
// r -> Random IC
// n -> Gaussian Normal IC
// u -> Uniform IC
char initCond = 'r';



// Level_set
class LS : public CF_3
{
public:
    LS() //this is the Level Set Function.
    {
        lip = 1.; //Lipschitz constant, for adaptive meshing purposes
    }
    
    double beta = 1.7; //radius of mother and daughter depend on this
    
    
    // ellipsoid constants //
    double a1 = 1.1;
    double a2 = 1.;
    
    double b1 = 1.2;
    double b2 = 1.;
    
    double c1 = 0.9;
    double c2 = 0.9;
    
    double x1_ell = xL/2.2; //where we center sphere
    double y1_ell = yL/2.2;
    double z1_ell = zL/2.2;
    
    double x2_ell = xL/2.2;
    double y2_ell = yL/2.2;
    double z2_ell = zL/2.2;
    
    // // //
    
    
    double r1_start = 0.15; //radius
    
    // so that we have mass conservation:
    double r2_start = 0.02; // + (0.005 * t); //another radius
    
    double alpha2 = (r1_start * (0.873580464) - r2_start)/9.0;
    
    double x1_start = xL/2.; //where we center sphere
    double y1_start = 0.25;
    double z1_start = zL/2.;
    
    double x2_start = xL/2.;
    
    double z2_start = zL/2.;
    
    
    tuple<double, double, double, double, double, double > Assign_phi(double x, double y, double z) const
    {
        
        double phi1, phi2, phi3, phi4, phi5, phi6;
        
        // (x-x1_start) *(y-y1_start) * (z-z1_start); //mother
        //double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t)); //daughter // sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*t);
        
        // Changed it to y1 start so that we start at the center of the mother cell
        // to make the daughter cell not exit the domain
        
        
        phi1 = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start);
        
        if (t < threshold)
        {
            phi2= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t));
        }
        else
        {
            phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*threshold)))+SQR(z-z2_start))-(r2_start+alpha2*(threshold));
            
        }
        //        double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start))+SQR(z-z2_start))-(r2_start);
        
        if (hole=='F')
        {
            return {phi1, phi2, phi3, phi4, phi5, phi6};
        }
        
        else if (hole=='T')
        {
            //this is where the level set comes in, where we're actually doing the level set:
            //we don't have negatives, so whichever gives us zero is our level set (because o(x) = 0, from picture)
            //double phi3 = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start/2);
            
            // phi3: Nucleus for the mother cell
            // phi4: Nucleus for the daughter cell
            
            
            if (t < 1.0)
            {
                phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start))+SQR(z-z1_start))-(r1_start/3);
                phi4= phi3/2;  //sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-((r2_start/3)+alpha2*(t));
            }
            
            else if (t >= 1.0 && t < 3.2)
            {
                phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start*(1+0.15*(t-1.0))))+SQR(z-z1_start))-(r1_start/3);
                phi4= phi3/2; //sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-((r2_start/3)+alpha2*(t));
            }
            
            else if (t > 3.2 && t < 7.6)
            {
                phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start*(1+0.15*(3.2-1.0))))+SQR(z-z1_start))-(r1_start/3);
                phi4= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*(t-0.9))))+SQR(z-z2_start))-((r2_start/4)+alpha2*(t/2.6));
            }
            
            else
            {
                phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start*(1+0.15*(9.8 - t))))+SQR(z-z1_start))-(r1_start/3);
                phi4= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-((r2_start/4)+alpha2*(t/2.6));
            }
            //            return MAX(MIN(phi1,phi2),-1*phi3);
            
            // when we want to enable the compartments
            if (compart == 'T')
            {
                phi5 = sqrt(SQR((x-x1_ell)/a1)+SQR((y-y1_ell)/b1)+SQR((z-z1_ell)/c1))-0.03;
                
                
                // this one we want to use or something that stays inside
                phi6 = sqrt(SQR((x-(x1_ell-0.01))/a1)+SQR((y-(y1_ell-0.17))/b1)+SQR((z-z1_ell)/c1))-0.03;
                
                // this one we want it to be outside
                //                phi6 = sqrt(SQR((x-(0.50))/a1)+SQR((y-(y1_ell))/b1)+SQR((z-z1_ell)/c1))-0.03;
                
                //              cout<< MIN(phi3,phi6) << "\n";
                return {phi1, phi2, phi3, phi4, phi5, phi6};;
                
            }
            
            // when we do not want additional compartments;=
            else if (compart == 'F')
            {
                return {phi1, phi2, phi3, phi4, phi5, phi6};;
            }
            
        }
        
        else
        {
            cout<<"Undefined Case of object in the middle...abort" << endl;
            return {phi1, phi2, phi3, phi4, phi5, phi6};
        }
        
    }
    
    
    double operator() (double x, double y, double z) const
    {
        
        auto [phi1, phi2, phi3, phi4, phi5, phi6] =  Assign_phi(x, y,z);
        
        
        if (hole=='F')
        {
            return MIN(phi1,phi2);
        }
        
        else if (hole=='T')
        {
            
            if (compart == 'T')
            {
                
                return MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi5,-1*phi6);
                
            }
            
            else if (compart == 'F')
            {
                return MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4));
            }
            
        }
        
        else
        {
            cout<<"Undefined Case of object in the middle...abort" << endl;
            return -1;
        }
        
    }
    
    
    int Find_Intersect(double x, double y, double z, int counter ,char whichOne)
    {
        auto [phi1, phi2, phi3, phi4, phi5, phi6] = Assign_phi(x,y,z);
        
        if (whichOne == '5')
        {
            // define all the geometry minus the desired level-set function
            
            
            double phi_total = MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi6);
            //
            if (sgn(phi5) == 1 && sgn(phi_total) == -1)
            {
                // means that they are colliding at some point
                counter++;
                
            }
            else
            {
                // do nothing for now
            }
            
            return counter;
        }
        
        // define all the geometry minus the desired level-set function
        else if (whichOne == '6')
        {
            double phi_total = MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi5);
            
            if (sgn(phi6) == -1 && sgn(phi_total) == 1)
            {
                // this compartment is colliding at some points
                counter++;
            }
            else
            {
                // do nothing for now
            }
            
            return counter;
        }
        
        else
        {
            // the right option was not provided
            return -1;
        }
    }
    
} level_set;



class Daughter_LvlSet : public CF_3
{
public:
    Daughter_LvlSet() //this is the Level Set Function.
    {
        lip = 1.; //Lipschitz constant, not important for us
    }
    double operator() (double x, double y, double z) const
    {
        
        double beta = 1.7; //radius of mother and daughter depend on this
        
        
        double r1_start = 0.15; //radius
        
        // so that we have mass conservation:
        double r2_start = 0.0; // + (0.005 * t); //another radius
        
        double alpha1 = r1_start/T*(beta-1); //placeholders
        double alpha2 = (r1_start * (0.873580464) - r2_start)/9.0;
        
        double x1_start = xL/2.; //where we center sphere
        double y1_start = 0.25;
        double z1_start = zL/2.;
        
        double x2_start = xL/2.;
        double y2_start = y1_start+r1_start;
        //        double y2_start = 0.75;
        
        double z2_start = zL/2.;
        
        double phi1 = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start); // (x-x1_start) *(y-y1_start) * (z-z1_start); //mother
        //double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t)); //daughter // sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*t);
        
        // Changed it to y1 start so that we start at the center of the mother cell
        
        // to make the daughter cell not exit the domain
        double phi2;
        double phi4;
        if (t < threshold)
        {
            phi2= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t));
        }
        else
        {
            phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*threshold)))+SQR(z-z2_start))-(r2_start+alpha2*(threshold));
            
        }
        //        double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start))+SQR(z-z2_start))-(r2_start);
        
        if (hole=='F')
        {
            return phi2;
        }
        
        else if (hole=='T')
        {
            //this is where the level set comes in, where we're actually doing the level set:
            //we don't have negatives, so whichever gives us zero is our level set (because o(x) = 0, from picture)
            
            // add a moving nucleous to the cells
            //            phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start*(1-0.15*(t-3.8))))+SQR(z-z1_start))-(r1_start/3);
            if (t > 4.8)
            {
                phi4= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-((r2_start/4)+alpha2*(t/2));
                return MAX(phi2,-1*phi4);
                
            }
            
        }
        // does nothing but just in case accidentally havea  different char for hole than F or T
        else
        {
            return phi2;
        }
        
        //this is where the level set comes in, where we're actually doing the level set: we don't have negatives, so whichever gives us zero is our level set (because o(x) = 0, from picture)
    }
} level_set_daughter;


class SplitCriteria : public SplitCriteriaOcTree //split criteria for the adaptive meshing
{
public:
    
    bool operator ()(const OcTree& tr, CaslInt c) const
    {
        double x = tr.get_Cell(c).icenter();
        double y = tr.get_Cell(c).jcenter();
        double z = tr.get_Cell(c).kcenter();
        
        if (level_set(x,y,z) < 0.001)
            return true;
        return false;
        
    }
    
};
// these are defined in the cpp files Phi_Zeta_Init_Class.cpp
Initial_Solution psi_initial;
zeta_Initial_Solution zeta_initial;
// Neumann Boundary condition
class WallBcPsiType : public WallBC3D
{
public:
    BoundaryConditionType operator()( double x, double y , double z) const
    {
        return NEUMANN; //just calling the class from Maxime's solver (for next parts as well)
    }
} wall_psi_neumann_type;
class WallBcPsiValues : public CF_3
{
public:
    double operator()(double x, double y, double z) const
    {
        return 0.;
    }
} wall_psi_neumann_value;
class IntBcPsiType : public WallBC3D
{
public:
    BoundaryConditionType operator()( double x, double y , double z) const
    {
        return NEUMANN;
    }
} int_psi_neumann_type;
class IntBcPsiValues : public CF_3
{
public:
    double operator()(double x, double y, double z) const
    {
        return 0.;
    }
} int_psi_neumann_value;


// ZETA BOUNDARY CONDITION
class WallBcZetaType : public WallBC3D
{
public:
    BoundaryConditionType operator()( double x, double y, double z) const
    {
        return NEUMANN;
    }
} wall_zeta_neumann_type;
class WallBcZetaValues : public CF_3
{
public:
    double operator()(double x, double y, double z) const
    {
        return 0.;
    }
} wall_zeta_neumann_value;
class IntBcZetaType : public WallBC3D
{
public:
    BoundaryConditionType operator()(double x, double y, double z) const
    {
        return NEUMANN;
    }
} int_zeta_neumann_type;
class IntBcZetaValues : public CF_3
{
    double operator()(double x, double y, double z) const
    {
        return 0.;
    }
} int_zeta_neumann_value;
//everything we had before was just declaring classes, but we now need to initialize them

int main(int argc, char **argv)
{
    /* PARALLALIZATION WITH OMP */
#ifdef CASL_OPENMP
    (void) omp_set_dynamic(false);
    //    if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
    (void) omp_set_num_threads(4);
    cout << "OpenMP - Max number of threads : " << omp_get_max_threads() << endl;
#endif
    
#ifdef CASL_PETSC
    Petsc::init(argc, argv);
#endif
    
    /* ************************************ */
    // this are defined in the cpp files Make_Folders.cpp
    /* CREATE A DIRECTORY FOR THE OUTPUT BASED ON THE TASK */
    MakeFolder(hole,DiffRXN,DiffOnly,D_psi, D_zeta, gamma_AtoB, gamma_BtoA, max_level, min_level, tOrder,
               &txtPathA,&txtPathB,&daughterPathA,&daughterPathB,&motherPathB,&motherPathA, FolderPath,initCond,FullPath,&Return_FolderPAth);
    
    FullPath = Return_FolderPAth;
    
    
    char whichOne1= '5';
    char whichOne2= '6';
    
    int coll_counter1 , coll_counter2;
    
    
    
    
    
    /* INITIALIZE OCTREE AND INITIAL CONDITIONS */
    
    my_octree.set_Grid(0.,xL,0.,yL,0.,zL); //setting the grid
    my_octree.construct_Octree_From_Level_Function(level_set, min_level, max_level);
    my_octree.impose_Uniform_Grid_On_Interface(level_set,max_level,5);
    //    my_octree.construct_Uniform_Octree(max_level);
    my_octree.initialize_Neighbors(); //talking about numerical cells as neighbors of each other
    
    
    ArrayV <double>  level_set_n (my_octree.number_Of_Nodes());
    
    // THIS IS WHAT I ADDED TO DO INTEGRATION OVER THE DOMAIN
    ArrayV<double> daughter_cell(my_octree.number_Of_Nodes());
    
    level_set_n.resize_Without_Copy(my_octree.number_Of_Nodes());
    daughter_cell.resize_Without_Copy(my_octree.number_Of_Nodes());
    
    psi_n.resize_Without_Copy(my_octree.number_Of_Leaves());
    psi_nm1.resize_Without_Copy(my_octree.number_Of_Leaves());
    
    
    //ZETA
    zeta_n.resize_Without_Copy(my_octree.number_Of_Leaves());
    zeta_nm1.resize_Without_Copy(my_octree.number_Of_Leaves());
    
    
    /************************************************************/
    // GENERATE THE INITIAL CONDITION HERE!!!
    
    IC_Generator init_cond;
    
    init_cond = gen(initCond,psi_n,zeta_n, my_octree,psi_initial,zeta_initial);
    psi_n = init_cond.psi;
    zeta_n = init_cond.zeta;
    
    /************************************************************/
    
#pragma omp parallel for
    for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
    {
        double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
        double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
        double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
        level_set_n(i) = level_set(x,y,z);
        
        // find the the two compartments are touching or nah
        coll_counter1 = level_set.Find_Intersect(x,y,z,coll_counter1,whichOne1);
        coll_counter2 = level_set.Find_Intersect(x,y,z,coll_counter2,whichOne2);
        
        daughter_cell(i) = level_set_daughter(x,y,z);
    }
    // to reinitialize the level-set and to take care of borderline cases
    OcTreeLevelSet ls_nodes(my_octree,level_set_n);
    ls_nodes.reinitialize();
    ls_nodes.perturb_Level_Function(1E-2*my_octree.dx_finest_resolution());
    ls_nodes.set_Phi(daughter_cell);
    ls_nodes.reinitialize();
    ls_nodes.perturb_Level_Function(1E-2*my_octree.dx_finest_resolution());
    
    
    
    /***************************************************/
    // print masses to files... here we initialize each file
    cout<<"PATH TXT: "<< txtPathA.c_str() <<endl;
    ofstream filestream(txtPathA.c_str());
    ofstream filestream2(txtPathB.c_str());
    ofstream daughterStreamA(daughterPathA.c_str());
    ofstream daughterStreamB(daughterPathB.c_str());
    ofstream motherStreamA(motherPathA.c_str());
    ofstream motherStreamB(motherPathB.c_str());
    
    /***************************************************/
    
    // MAIN WHILE LOOP
    
    while (t <= T)
    {
        
        cout<<"We are at step: "<<n<<" at time = "<<t<<endl;
        OcTreeSolverCellCenteredPoisson poisson_solver;
        poisson_solver.set_Octree(my_octree);
        BoundaryConditions3D bc_psi;
        bc_psi.setWallTypes(wall_psi_neumann_type);
        bc_psi.setWallValues(wall_psi_neumann_value);
        bc_psi.setInterfaceType(NEUMANN);
        bc_psi.setInterfaceValue(int_psi_neumann_value);
        
        
        cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << coll_counter1 << " many collsions of Phi " << whichOne1 <<" compartments detected" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << coll_counter2 << " many collsions of Phi " << whichOne2 <<" compartments detected" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
        cout<<"Resetting counters to 0 " << endl;
        // reset counters
        coll_counter1 = 0;
        coll_counter2 = 0;
        
        
        
        
        // at time step 0 we do not want to solve anything, just intital conditions
        if (n > 0)
        {
            ArrayV<double> psi_updated;
            psi_updated.resize_Without_Copy(my_octree.number_Of_Leaves());
            
            //#pragma omp parallel for
            
            
            
            for (int i = 0; i < psi_updated.size(); i++) //because we got adding array operator, just inbetween step
            {
                // rhs
                // for now we are doing Diffusion Only
                psi_updated(i) =psi_n(i) + (initial_pop*dt) - (dt*gamma_AtoB*(psi_n(i))*(psi_n(i))) + (2 *dt* gamma_BtoA * zeta_n(i));// MOST RECENT -> DIFFUSION -> psi_n(i) + 2*gamma_BtoA * zeta_nm1(i); // (alpha - beta*(psi_n(i))*(psi_n(i)) + gamma_BtoA * zeta_m(i));//2 * gamma_BtoA * zeta_m(i))
                
                
            }
            
            
            ArrayV<double> rhs = psi_updated;//previous right hand side, kind of like a placeholder
            // here we will use the previous volume to integrate the rhs for mass conservation
            ArrayV<double> level_set_nm1(my_octree.number_Of_Nodes());
            t-=dt;
#pragma omp parallel for
            for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
            {
                double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
                double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
                double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
                level_set_nm1(i) = level_set(x,y,z);
            }
            ls_nodes.set_Phi(level_set_nm1);
            ls_nodes.set_Octree(my_octree);
            ls_nodes.reinitialize();
            ls_nodes.perturb_Level_Function(1E-2*my_octree.dx_finest_resolution());
            t+=dt;
            rhs.CHK_NAN();
#pragma omp parallel for
            for (int i = 0; i < my_octree.number_Of_Leaves(); i++) //because we got adding array operator, just inbetween step
            {
                CaslInt c  = my_octree.leaf2cell(i);
                const OctCell& C = my_octree.get_Cell(c);
                // need to find the Vin and Vinp1 volumes
                Cube3 cube;
                // this is the cell C_i
                cube.x0 = my_octree.x_fr_i(C.imin());
                cube.x1 = my_octree.x_fr_i(C.imax());
                cube.y0 = my_octree.y_fr_j(C.jmin());
                cube.y1 = my_octree.y_fr_j(C.jmax());
                cube.z0 = my_octree.z_fr_k(C.kmin());
                cube.z1 = my_octree.z_fr_k(C.kmax());
                
                // this is the value of the level set a the corners of the cube
                OctValue leveset_nm1_values(level_set_nm1(my_octree.get_Cell(c).node_mmm()),level_set_nm1(my_octree.get_Cell(c).node_mmp()),
                                            level_set_nm1(my_octree.get_Cell(c).node_mpm()),level_set_nm1(my_octree.get_Cell(c).node_mpp()),
                                            level_set_nm1(my_octree.get_Cell(c).node_pmm()),level_set_nm1(my_octree.get_Cell(c).node_pmp()),
                                            level_set_nm1(my_octree.get_Cell(c).node_ppm()),level_set_nm1(my_octree.get_Cell(c).node_ppp()));
                
                double Vnm1 = cube.volume_In_Negative_Domain(leveset_nm1_values);
                OctValue leveset_n_values(level_set_n(my_octree.get_Cell(c).node_mmm()),level_set_n(my_octree.get_Cell(c).node_mmp()),
                                          level_set_n(my_octree.get_Cell(c).node_mpm()),level_set_n(my_octree.get_Cell(c).node_mpp()),
                                          level_set_n(my_octree.get_Cell(c).node_pmm()),level_set_n(my_octree.get_Cell(c).node_pmp()),
                                          level_set_n(my_octree.get_Cell(c).node_ppm()),level_set_n(my_octree.get_Cell(c).node_ppp()));
                double Vn = cube.volume_In_Negative_Domain(leveset_n_values);
                
                //integrate and divide by future volume Vn so that its cancel after solver integrate over Vn
                
                rhs(i) *=Vnm1/MAX(EPSILON,Vn);
                
                
            }
            
            rhs.CHK_NAN();
            cout<<rhs.max_Abs()<<endl;
            cout<<rhs.min()<<endl;
            //ArrayV<double> rhs = psi_nm1;
            poisson_solver.set_bc(bc_psi);
            poisson_solver.set_Rhs(rhs);
            // the diagonals here are what need to be changed if we want to add reaction terms
            poisson_solver.set_Diagonal_Increment(1.0 + mu * dt);
            
            poisson_solver.set_mu(dt*D_psi);
            
            
            poisson_solver.set_Phi(level_set_n);
            poisson_solver.set_Linear_Solver(solver_PETSC);
            poisson_solver.solve(psi_n); //we then solve for psi_n
            psi_n.CHK_NAN();
        }
        
        
        
        OcTreeCellBasedLevelSet ls;
        ls.set_Octree(my_octree);
        ls.set_Phi(level_set_n);
        
        OcTreeCellBasedLevelSet ls_daughter;
        ls_daughter.set_Octree(my_octree);
        ls_daughter.set_Phi(daughter_cell);
        
        
        double total_psi = ls.integral_Cell_Based(psi_n);
        double total_zeta = ls.integral_Cell_Based(zeta_n);
        
        double only_daughter_psi = ls_daughter.integral_Cell_Based(psi_n);
        double only_daughter_zeta = ls_daughter.integral_Cell_Based(zeta_n);
        
        
        
        if (n==0)
        {
            init_mass_A = total_psi;
            // since we have A + A = B
            init_mass_B = 2 * total_zeta;
            
        }
        
        
        //ZETA
        BoundaryConditions3D bc_zeta;
        bc_zeta.setWallTypes(wall_zeta_neumann_type);
        bc_zeta.setWallValues(wall_zeta_neumann_value);
        bc_zeta.setInterfaceType(NEUMANN);
        bc_zeta.setInterfaceValue(int_zeta_neumann_value);
        
        // at time step 0 we do not want to solve anything, just intital conditions
        if (n > 0)
        {
            
            ArrayV<double> zeta_updated;
            zeta_updated.resize_Without_Copy(my_octree.number_Of_Leaves());
            
            //#pragma omp parallel for
            
            
            for (int i = 0; i < zeta_updated.size(); i++) //because we got adding array operator, just inbetween step
            {
                // rhs
                // since zeta is a dimer, then two monomers give you one dimer so that's why there's the 0.5 in there
                zeta_updated(i) = zeta_n(i) + 0.5 * (dt*gamma_AtoB*(psi_n(i))*(psi_n(i)));
                // included in the LHS of the solver - (dt * gamma_BtoA * zeta_n(i));
                
                
            }
            
            
            ArrayV<double> rhs_z =zeta_updated;
            ArrayV<double> level_set_nm1(my_octree.number_Of_Nodes());
            t-=dt;
#pragma omp parallel for
            for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
            {
                double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
                double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
                double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
                level_set_nm1(i) = level_set(x,y,z);
            }
            ls_nodes.set_Phi(level_set_nm1);
            ls_nodes.set_Octree(my_octree);
            ls_nodes.reinitialize();
            ls_nodes.perturb_Level_Function(1E-2*my_octree.dx_finest_resolution());
            t+=dt;
            
#pragma omp parallel for
            for (int i = 0; i < my_octree.number_Of_Leaves(); i++) //because we got adding array operator, just inbetween step
            {
                CaslInt c  = my_octree.leaf2cell(i);
                const OctCell& C = my_octree.get_Cell(c);
                // need to find the Vin and Vinp1 volumes
                Cube3 cube;
                // this is the cell C_i
                cube.x0 = my_octree.x_fr_i(C.imin());
                cube.x1 = my_octree.x_fr_i(C.imax());
                cube.y0 = my_octree.y_fr_j(C.jmin());
                cube.y1 = my_octree.y_fr_j(C.jmax());
                cube.z0 = my_octree.z_fr_k(C.kmin());
                cube.z1 = my_octree.z_fr_k(C.kmax());
                
                // this is the value of the level set a the corners of the cube
                OctValue leveset_nm1_values(level_set_nm1(my_octree.get_Cell(c).node_mmm()),level_set_nm1(my_octree.get_Cell(c).node_mmp()),
                                            level_set_nm1(my_octree.get_Cell(c).node_mpm()),level_set_nm1(my_octree.get_Cell(c).node_mpp()),
                                            level_set_nm1(my_octree.get_Cell(c).node_pmm()),level_set_nm1(my_octree.get_Cell(c).node_pmp()),
                                            level_set_nm1(my_octree.get_Cell(c).node_ppm()),level_set_nm1(my_octree.get_Cell(c).node_ppp()));
                
                double Vnm1 = cube.volume_In_Negative_Domain(leveset_nm1_values);
                OctValue leveset_n_values(level_set_n(my_octree.get_Cell(c).node_mmm()),level_set_n(my_octree.get_Cell(c).node_mmp()),
                                          level_set_n(my_octree.get_Cell(c).node_mpm()),level_set_n(my_octree.get_Cell(c).node_mpp()),
                                          level_set_n(my_octree.get_Cell(c).node_pmm()),level_set_n(my_octree.get_Cell(c).node_pmp()),
                                          level_set_n(my_octree.get_Cell(c).node_ppm()),level_set_n(my_octree.get_Cell(c).node_ppp()));
                double Vn = cube.volume_In_Negative_Domain(leveset_n_values);
                
                //integrate and divide by future volume Vn so that its cancel after solver integrate over Vn
                
                rhs_z(i) *=Vnm1/MAX(EPSILON,Vn);
                
                
            }
            poisson_solver.set_bc(bc_zeta);
            poisson_solver.set_Rhs(rhs_z);
            poisson_solver.set_Linear_Solver(solver_PETSC);
            // the diagonals here are what need to be changed if we want to add reaction terms
            poisson_solver.set_Diagonal_Increment(1.0 + (mu * dt) + (dt * gamma_BtoA));
            poisson_solver.set_mu(dt*D_zeta);
            
            poisson_solver.set_Phi(level_set_n);
            poisson_solver.solve(zeta_n);
            
        }
        //        zeta_n.CHK_NAN();
        //        psi_n.CHK_NAN();
        double mother_psi = total_psi - only_daughter_psi ;
        double mother_zeta =  total_zeta - only_daughter_zeta ;
        
        
        cout << "Psi in Daughter: " << only_daughter_psi << endl;
        cout << "Zeta in Daughter: " << only_daughter_zeta << endl;
        
        cout << "total mass in daughter cell : "<< only_daughter_psi + 2* only_daughter_zeta << endl;
        cout << "total mass in mother cell : "<< mother_psi + 2 * mother_zeta << endl;
        
        // No need for this with the our new FV method
        //            ls.extrapolate_Along_Normal_Using_Cells(psi_n, bc_psi);
        //            ls.extrapolate_Along_Normal_Using_Cells(zeta_n, bc_zeta);
        
        total_psi = ls.integral_Cell_Based(psi_n);
        total_zeta = ls.integral_Cell_Based(zeta_n);
        
        // function call for detecting NaNs and throwing exceptions
        /////////////////////////////
        IsNan(total_psi,total_zeta);
        /////////////////////////////
        
        cout << "Current Total Psi: "<< total_psi << endl;
        cout << "Current Total zeta: "<< total_zeta << endl;
        
        // individual masses of species
        filestream << n << "     " << total_psi << endl;
        filestream2 << n << "     " << total_zeta << endl;
        
        // individual masses of species in daughter cell
        daughterStreamA << n << "     " << only_daughter_psi << endl;
        daughterStreamB<< n << "     " << only_daughter_zeta << endl;
        
        // individual masses of species in mother cell
        motherStreamA << n << "     " << mother_psi << endl;
        motherStreamB<< n << "     " << mother_zeta << endl;
        
        // we have A + A = B -> B = 2A
        cout << "***********************************" << endl;
        cout << "Current Total mass: " << total_psi + 2 * total_zeta << endl;
        cout << "-----------------------------------" << endl;
        
        t += dt;
        
        
        
        new_octree.make_It_Have_Root_Cell_Only();
        new_octree.set_Grid(0.,xL,0.,yL,0.,zL);
        new_octree.construct_Octree_From_Level_Function(level_set, min_level, max_level);
        new_octree.impose_Uniform_Grid_On_Interface(level_set,max_level,5);
        
        new_octree.initialize_Neighbors();
        
        
        
        // update the quatity in the new mesh
        ArrayV <double> new_psi_n (new_octree.number_Of_Leaves());
        //ZETA
        ArrayV <double> new_zeta_n(new_octree.number_Of_Leaves());
#pragma omp parallel for
        for (CaslInt leaf = 0; leaf<new_octree.number_Of_Leaves(); leaf++)
        {
            double x = new_octree.x_fr_i(new_octree.get_Leaf(leaf).icenter());
            double y = new_octree.y_fr_j(new_octree.get_Leaf(leaf).jcenter());
            double z = new_octree.z_fr_k(new_octree.get_Leaf(leaf).kcenter());
            
            new_psi_n(leaf) = interpolation_LSQR_Cell_Centered(my_octree,psi_n, x, y, z, NEUMANN, level_set_n);
            //ZETA
            new_zeta_n(leaf) = interpolation_LSQR_Cell_Centered(my_octree, zeta_n, x, y, z, NEUMANN, level_set_n);
        }
        // for visu only
        my_octree = new_octree;
        psi_nm1 = new_psi_n;
        psi_n = new_psi_n;
        //ZETA
        zeta_nm1 = new_zeta_n;
        zeta_n = new_zeta_n;
        // to check nans
        //        psi_n.CHK_NAN();
        //        zeta_n.CHK_NAN();
        
        
        level_set_n.resize_Without_Copy(my_octree.number_Of_Nodes());
        daughter_cell.resize_Without_Copy(my_octree.number_Of_Nodes());
#pragma omp parallel for
        for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
        {
            double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
            double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
            double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
            level_set_n(i) = level_set(x,y,z);
            
            // check if the compartments are touching as we update the level-set
            coll_counter1 = level_set.Find_Intersect(x,y,z,coll_counter1,whichOne1);
            coll_counter2 = level_set.Find_Intersect(x,y,z,coll_counter1,whichOne2);
            daughter_cell(i) = level_set_daughter(x,y,z);
            
        }
        
        // to reinitialize the level-set and to take care of borderline cases
        ls_nodes.set_Phi(level_set_n);
        ls_nodes.set_Octree(my_octree);
        ls_nodes.reinitialize();
        double dx,dy,dz;
        my_octree.dx_dy_dz_smallest(dx,dy,dz);
        ls_nodes.perturb_Level_Function(1E-2*dx);
        ls_nodes.set_Phi(daughter_cell);
        ls_nodes.reinitialize();
        ls_nodes.perturb_Level_Function(1E-2*dx);
        ArrayV<double> psi_node(my_octree.number_Of_Nodes());
        //ZETA
        ArrayV<double> zeta_node(my_octree.number_Of_Nodes());
#pragma omp parallel for
        for (CaslInt n = 0; n<my_octree.number_Of_Nodes(); n++)
        {
            double x = my_octree.x_fr_i(my_octree.get_Node(n).i);
            double y = my_octree.y_fr_j(my_octree.get_Node(n).j);
            double z = my_octree.z_fr_k(my_octree.get_Node(n).k);
            
            psi_node(n) = interpolation_LSQR_Cell_Centered(my_octree, psi_n, x, y, z, NEUMANN, level_set_n);
            //ZETA
            zeta_node(n) = interpolation_LSQR_Cell_Centered(my_octree, zeta_n, x, y, x, NEUMANN, level_set_n);
        }
        
        
        char file_name [500];
        
        sprintf(file_name,  "%s/cell_splitting_%d.vtk",FullPath.c_str(),n);
        my_octree.print_VTK_Format(file_name);
        my_octree.print_VTK_Format(level_set_n, "level_set",file_name);
        my_octree.print_VTK_Format(psi_node, "psi_node",file_name);
        //ZETA
        my_octree.print_VTK_Format(zeta_node, "zeta_node",file_name);
        my_octree.print_VTK_Format_Switch_To_Cell(file_name);
        my_octree.print_VTK_Format_Cell(psi_n, "psi",file_name);
        //ZETA
        my_octree.print_VTK_Format_Cell(zeta_n, "zeta",file_name);
        
        
        n++;
        
    }
    
    return 0;
}
