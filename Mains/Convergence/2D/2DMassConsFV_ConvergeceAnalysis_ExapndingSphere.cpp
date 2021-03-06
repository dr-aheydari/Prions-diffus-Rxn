#ifdef CASL_OPENMP
#include <omp.h>
#define TRUE  1
#define FALSE 0
#else
#define omp_get_thread_num() 0
#endif


#include <list>
#include <tuple>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/uio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include </Users/aliheydari/stdc++.h>


// Maxime's proprietary libraries - NOT included in the repos for obvious reasons
#include <lib/amr/QuadNode.h>
#include <lib/amr/QuadCell.h>
#include <lib/amr/QuadNgbdNodesOfNode.h>

#include <lib/arrays/ArrayV.h>
#include <lib/amr/QuadTree.h>
#include <lib/amr/QuadTreeCellBased.h>
#include <lib/amr/QuadTreeLevelSet.h>
#include <lib/amr/QuadTreeAdvection.h>
#include <lib/tools/QuadTreeToolbox.h>
#include <lib/amr/QuadTreeCellBasedLevelSet.h>
#include <lib/solvers/QuadTreeSolverCellCenteredPoisson.h>

// My custom tools - included in the final project repository for reproducibility
#include <lib/AliTools/isNan.h>
#include <lib/AliTools/vectorUtils.h>
#include <lib/AliTools/IC_Generator.h>
#include <lib/AliTools/Make_Folders.h>

#include <lib/AliTools/signum.cpp>
#include <lib/AliTools/Psi_Zeta_Init_Class.cpp>

#ifdef CASL_PETSC
#include <lib/tools/PetscTools.h>
#endif


using namespace std;
using namespace CASL;

////////////
int min_level = 9; //for our adaptive meshing
int max_level = 13; //for our adaptive meshing
////////////

//gives a coordinate system
QuadTreeCellBased my_quadtree;
QuadTreeCellBased new_quadtree;

ArrayV <double> psi_n;
ArrayV <double> psi_nm1;

//ZETA VARIABLES
ArrayV <double> zeta_n;
ArrayV <double> zeta_nm1;

//// Space discretization of [0,L]
double xL = 0.75;
double yL = 0.75;
double zL = 0.75;


// Time discretization
double T = 5; //final time
double t = 0; // initial time
double dt = 0.1; //time step
double tOrder = dt;
int n = 0; // just for iterations
double dt_const = 1000.; // the fixed constant that multiplies dt for speed


// System parameters
double D_psi = 1; // Diffusion coeffcient for "healthy" protein
double D_zeta = 0.00; // Diffusion coeffcient for aggregate
//// for convergence study we only look at one species

//rates
double gamma_AtoB = 0.001 * 0;
double initial_pop = 10*0;
double gamma_BtoA = 0.01 * 0;
double mu = 0.2 * 0;

double init_mass_A;

double epsilon = 0.000000001;

// this threshold is used for knowing when the daughter cell will leave the mother cell
double threshold = 9;

double vol_n, err_t = 0;

// intialize for explicit task dependent file IO
char hole = 'F';
char DiffRXN = 'T';
char DiffOnly='F';
char compart = 'F';

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
// e -> exact IC for the exact analytical solution (done in convergence)
char initCond = 'e';

// Level_set
class LS : public CF_2
{
public:
    LS() //this is the Level Set Function.
    {
        lip = 1.; //Lipschitz constant, for adaptive meshing purposes
    }
    
    double beta = 1.7; //growth factor-radius of mother and daughter depend on this
    
    // // //
    
    double r1_start = 0.05; // radius for the changing domain
//    double r1_start = 0.25; //radius for fixed domain
    
    // so that we have mass conservation:
    double r2_start = 0.02; // + (0.005 * t); //another radius
    
    double alpha2 = (r1_start * (0.873580464) - r2_start)/2;
    
    double x1_start = xL/2.; //where we center sphere
    double y1_start = 0.25;
    
    double Assign_phi(double x, double y) const
    {
        
        double phi1;
        
        // changing
        phi1 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start*(1+0.05*t))))-(r1_start+alpha2*t);
        return phi1;
            
    }
    
    double operator() (double x, double y) const
    {

        double phi1 = Assign_phi(x,y);
        return phi1;

    }
    
} level_set;



class SplitCriteria : public SplitCriteriaQuadTree //split criteria for the adaptive meshing
{
public:
    
    bool operator ()(const QuadTree& tr, CaslInt c) const
    {
        double x = tr.get_Cell(c).icenter();
        double y = tr.get_Cell(c).jcenter();
        
        if (level_set(x,y) < 0.001)
            return true;
        return false;
        
    }
    
};

Initial_Solution2D psi_initial;
zeta_Initial_Solution2D zeta_initial;
// Neumann Boundary condition
class WallBcPsiType : public WallBC2D
{
public:
    BoundaryConditionType operator()( double x, double y) const
    {
        return NEUMANN; //just calling the class from Maxime's solver (for next parts as well)
    }
} wall_psi_neumann_type;
class WallBcPsiValues : public CF_2
{
public:
    double operator()(double x, double y) const
    {
        return 0;
    }

} wall_psi_neumann_value;
class IntBcPsiType : public WallBC2D
{
public:
    BoundaryConditionType operator()( double x, double y) const
    {
        return NEUMANN;
    }
} int_psi_neumann_type;

// interface boundary condition
class IntBcPsiValues : public CF_2
{
public:
    double operator()(double x, double y) const
    {
        /*
        double r1_start = 0.05; // radius for the changing domain
        double r2_start = 0.02; // + (0.005 * t); //another radius
        */

        double x1_start = xL/2.; //where we center sphere
        double y1_start = 0.25;
        // gradient of the level set
        double ls_grad = (x-x1_start) + (y-(y1_start*(1+0.05*t)));

        double r = sqrt(SQR(x-x1_start)+SQR(y-(y1_start*(1+0.05*t))));
        double psi_exact = exp(x + y + t);

        return (1/r) * psi_exact * ls_grad;
    }

} int_psi_neumann_value;




int main(int argc, char **argv)
{
    /* PARALLALIZATION WITH OMP */
#ifdef CASL_OPENMP
    (void) omp_set_dynamic(false);
    if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
    (void) omp_set_num_threads(4);
    cout << "OpenMP - Max number of threads : " << omp_get_max_threads() << endl;
#endif
    
#ifdef CASL_PETSC
    Petsc::init(argc, argv);
#endif
    
    /* ************************************ */
    // this are defined in the cpp files Make_Folders.cpp
    /* CREATE A DIRECTORY FOR THE OUTPUT BASED ON THE TASK */


    MakeFolder(compart, hole,DiffRXN,DiffOnly,D_psi, D_zeta, gamma_AtoB, gamma_BtoA, max_level, min_level, tOrder,
               &txtPathA,&txtPathB,&daughterPathA,&daughterPathB,&motherPathB,&motherPathA, FolderPath,initCond,FullPath,&Return_FolderPAth);
    
    FullPath = Return_FolderPAth;
    
    /* INITIALIZE OCTREE AND INITIAL CONDITIONS */
    
    my_quadtree.set_Grid(0.,xL,0.,yL); //setting the grid
    my_quadtree.construct_Quadtree_From_Level_Function(level_set, min_level, max_level);
    my_quadtree.impose_Uniform_Grid_On_Interface(level_set,max_level,5);

    my_quadtree.initialize_Neighbors(); //talking about numerical cells as neighbors of each other
    
    
    ArrayV <double>  level_set_n (my_quadtree.number_Of_Nodes());
    ArrayV <double> psi_exact (my_quadtree.number_Of_Nodes());
    
    level_set_n.resize_Without_Copy(my_quadtree.number_Of_Nodes());
    psi_exact.resize_Without_Copy(my_quadtree.number_Of_Nodes());
    
    psi_n.resize_Without_Copy(my_quadtree.number_Of_Leaves());
    psi_nm1.resize_Without_Copy(my_quadtree.number_Of_Leaves());
    
    
    
    /************************************************************/
    // GENERATE THE INITIAL CONDITION HERE!!!
    
    IC_Generator init_cond;
    init_cond = gen_exact2D(initCond,psi_n, my_quadtree,psi_initial);
    psi_n = init_cond.psi;
    
    /************************************************************/
    
#pragma omp parallel for
    for (int i = 0; i<my_quadtree.number_Of_Nodes(); i++)
    {
        double x = my_quadtree.x_fr_i(my_quadtree.get_Node(i).i);
        double y = my_quadtree.y_fr_j(my_quadtree.get_Node(i).j);
        level_set_n(i) = level_set(x,y);
        psi_exact(i) = exp(x + y + t);

    }

    psi_n = psi_exact;
    // to reinitialize the level-set and to take care of borderline cases
    QuadTreeLevelSet ls_nodes(my_quadtree,level_set_n);
    ls_nodes.reinitialize();
    ls_nodes.perturb_Level_Function(1E-2*my_quadtree.dx_finest_resolution());

    ls_nodes.set_Phi(level_set_n);
    ls_nodes.reinitialize();
    ls_nodes.perturb_Level_Function(1E-2*my_quadtree.dx_finest_resolution());
    
    
    
    /***************************************************/
    // print masses to files... here we initialize each file
    cout<<"PATH TXT: "<< txtPathA.c_str() <<endl;
    ofstream filestream(txtPathA.c_str());

    
    /***************************************************/
    
    // MAIN WHILE LOOP
    while (t <= T)
    {

        cout<<"Here is my dt: " << dt << endl;
        cout<<"We are at step: "<<n<<" at time = "<<t<<endl;
        QuadTreeSolverCellCenteredPoisson poisson_solver;
        poisson_solver.set_Quadtree(my_quadtree);
        BoundaryConditions2D bc_psi;
        bc_psi.setWallTypes(wall_psi_neumann_type);
        bc_psi.setWallValues(wall_psi_neumann_value);
        bc_psi.setInterfaceType(NEUMANN);
        bc_psi.setInterfaceValue(int_psi_neumann_value);   
        
        double total_psi;
        double total_psi_exact;

        QuadTreeCellBasedLevelSet ls;
        ls.set_Quadtree(my_quadtree);
        ls.set_Phi(level_set_n);

        if (n==0)
            {
                total_psi = ls.integral_Cell_Based(psi_n);
                total_psi_exact = ls.integral_Cell_Based(psi_exact);

                ArrayV<double> psi_node(my_quadtree.number_Of_Nodes());
                ArrayV<double> psi_exact_node(my_quadtree.number_Of_Nodes());

    #pragma omp parallel for
            for (CaslInt n = 0; n<my_quadtree.number_Of_Nodes(); n++)
                {
                    double x = my_quadtree.x_fr_i(my_quadtree.get_Node(n).i);
                    double y = my_quadtree.y_fr_j(my_quadtree.get_Node(n).j);

                    psi_node(n) = interpolation_LSQR_Cell_Centered(my_quadtree, psi_n, x, y, NEUMANN, level_set_n);
                    psi_exact_node(n) = interpolation_LSQR_Cell_Centered(my_quadtree, psi_exact, x, y, NEUMANN, level_set_n);
                }


                cout << "***********************************" << endl;
                cout << "Current Total mass: " << total_psi << endl;
                cout << "-----------------------------------" << endl;

                char file_name [500];

                sprintf(file_name,  "%s/cell_splitting_%d.vtk",FullPath.c_str(),n);
                my_quadtree.print_VTK_Format(file_name);
                my_quadtree.print_VTK_Format(level_set_n, "level_set",file_name);
                my_quadtree.print_VTK_Format(psi_node, "psi_node",file_name);
                my_quadtree.print_VTK_Format(psi_exact_node, "psi_exact_node",file_name);
                t += dt;

                cout<< "+++++++++++++++++++++++++" << endl;
                cout<<"VOLUME: "<< ls_nodes.area_In_Negative_Domain() << endl;
                err_t = ABS(total_psi - total_psi_exact)/pow(2,max_level);
                cout<<"Error:"<< err_t <<endl;
                cout<< "+++++++++++++++++++++++++" << endl;

            }

        // at time step 0 we do not want to solve anything, just intital conditions
        if (n > 0)
            {
                ArrayV<double> psi_updated;
                psi_updated.resize_Without_Copy(my_quadtree.number_Of_Leaves());

        #pragma omp parallel for
                for (int i = 0; i < psi_updated.size(); i++) //because we got adding array operator, just inbetween step
                {
                    // rhs
                    psi_updated(i) = psi_n(i) + (dt * -2 * psi_n(i));
                }

                ArrayV<double> rhs = psi_updated;
                // here we will use the previous volume to integrate the rhs for mass conservation
                ArrayV<double> level_set_nm1(my_quadtree.number_Of_Nodes());

                t-=dt;
        #pragma omp parallel for
                for (int i = 0; i<my_quadtree.number_Of_Nodes(); i++)
                    {
                        double x = my_quadtree.x_fr_i(my_quadtree.get_Node(i).i);
                        double y = my_quadtree.y_fr_j(my_quadtree.get_Node(i).j);
                        level_set_nm1(i) = level_set(x,y);
                    }
                ls_nodes.set_Phi(level_set_nm1);
                ls_nodes.set_Quadtree(my_quadtree);
                ls_nodes.reinitialize();
                ls_nodes.perturb_Level_Function(1E-2*my_quadtree.dx_finest_resolution());

                t+=dt;
                rhs.CHK_NAN();

/* This is our modification to FV in 2D
 
#pragma omp parallel for
                for (int i = 0; i < my_quadtree.number_Of_Leaves(); i++) //because we got adding array operator, just inbetween step
                {
                    CaslInt c  = my_quadtree.leaf2cell(i);
                    const QuadCell& C = my_quadtree.get_Cell(c);
                    
                    // need to find the Vin and Vinp1 volumes
                    Cube2 cube;
                    // this is the cell C_i
                    cube.x0 = my_quadtree.x_fr_i(C.imin());
                    cube.x1 = my_quadtree.x_fr_i(C.imax());
                    cube.y0 = my_quadtree.y_fr_j(C.jmin());
                    cube.y1 = my_quadtree.y_fr_j(C.jmax());
                    
                    // this is the value of the level set a the corners of the cube
                    QuadValue leveset_nm1_values(level_set_nm1(my_quadtree.get_Cell(c).node_mm()),level_set_nm1(my_quadtree.get_Cell(c).node_mp()),
                                                 level_set_nm1(my_quadtree.get_Cell(c).node_pm()),level_set_nm1(my_quadtree.get_Cell(c).node_pp())
                                                 );
                    
                    double Vnm1 = cube.area_In_Negative_Domain(leveset_nm1_values);
                    
                    QuadValue leveset_n_values(level_set_n(my_quadtree.get_Cell(c).node_mm()),level_set_n(my_quadtree.get_Cell(c).node_mp()),
                                               level_set_n(my_quadtree.get_Cell(c).node_pm()),level_set_n(my_quadtree.get_Cell(c).node_pp())
                                               );
                    
                    double Vn = cube.area_In_Negative_Domain(leveset_n_values);
                    //integrate and divide by future volume Vn so that its cancel after solver integrate over Vn
                    rhs(i) *=Vnm1/MAX(EPSILON,Vn);
                    
                }
*/
                

//////////////////////////////////////////////////// Debugging: check for any NaNs////////////////////////////////////////////////////////
//                rhs.CHK_NAN();
//                cout<<rhs.max_Abs()<<endl;
//                cout<<rhs.min()<<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                poisson_solver.set_bc(bc_psi);
                poisson_solver.set_Rhs(rhs);
                // the diagonals here are what need to be changed if we want to add reaction terms
                poisson_solver.set_Diagonal_Increment(1.0 + (mu * dt) + (dt * gamma_AtoB));

                poisson_solver.set_mu(dt*D_psi);
                poisson_solver.set_Phi(level_set_n);
                poisson_solver.set_Linear_Solver(solver_PETSC);
                poisson_solver.solve(psi_n); //we then solve for psi_n

                QuadTreeCellBasedLevelSet ls;
                ls.set_Quadtree(my_quadtree);
                ls.set_Phi(level_set_n);

                // to keep track of the total mass as function of time

                total_psi = ls.integral_Cell_Based(psi_n);
                total_psi_exact = ls.integral_Cell_Based(psi_exact);

                cout << "Current Total Psi: "<< total_psi << endl;
                cout << "Current Total Psi Exact: "<< total_psi_exact << endl;

                cout << "***********************************" << endl;
                cout << "Current Total mass: " << total_psi << endl;
                cout << "-----------------------------------" << endl;


                new_quadtree.make_It_Have_Root_Cell_Only();
                new_quadtree.set_Grid(0.,xL,0.,yL);
                new_quadtree.construct_Quadtree_From_Level_Function(level_set, min_level, max_level);
                new_quadtree.impose_Uniform_Grid_On_Interface(level_set,max_level,5);

                new_quadtree.initialize_Neighbors();

                // update the quatity in the new mesh
                ArrayV <double> new_psi_n (new_quadtree.number_Of_Leaves());
            #pragma omp parallel for
                for (CaslInt leaf = 0; leaf<new_quadtree.number_Of_Leaves(); leaf++)
                    {
                        double x = new_quadtree.x_fr_i(new_quadtree.get_Leaf(leaf).icenter());
                        double y = new_quadtree.y_fr_j(new_quadtree.get_Leaf(leaf).jcenter());

                        new_psi_n(leaf) = interpolation_LSQR_Cell_Centered(my_quadtree,psi_n, x, y, NEUMANN, level_set_n);
                    }

                // for visualization purposes
                my_quadtree = new_quadtree;
                psi_nm1 = new_psi_n;
                psi_n = new_psi_n;

                level_set_n.resize_Without_Copy(my_quadtree.number_Of_Nodes());
                psi_exact.resize_Without_Copy(my_quadtree.number_Of_Nodes());

            #pragma omp parallel for
                for (int i = 0; i<my_quadtree.number_Of_Nodes(); i++)
                    {
                        double x = my_quadtree.x_fr_i(my_quadtree.get_Node(i).i);
                        double y = my_quadtree.y_fr_j(my_quadtree.get_Node(i).j);
                        level_set_n(i) = level_set(x,y);
                        psi_exact(i) = exp(x+y+t);

                    }

                // to reinitialize the level-set and to take care of borderline cases
                ls_nodes.set_Phi(level_set_n);
                ls_nodes.set_Quadtree(my_quadtree);
                ls_nodes.reinitialize();
                double dx,dy;
                my_quadtree.dx_and_dy_smallest(dx,dy);
                ls_nodes.perturb_Level_Function(1E-2*dx);

                ArrayV<double> psi_node(my_quadtree.number_Of_Nodes());
                ArrayV<double> psi_exact_node(my_quadtree.number_Of_Nodes());
                ArrayV<double> sols_diff(my_quadtree.number_Of_Nodes());
                sols_diff.resize_Without_Copy(my_quadtree.number_Of_Nodes());



            #pragma omp parallel for
                for (CaslInt n = 0; n<my_quadtree.number_Of_Nodes(); n++)
                    {
                        double x = my_quadtree.x_fr_i(my_quadtree.get_Node(n).i);
                        double y = my_quadtree.y_fr_j(my_quadtree.get_Node(n).j);
                        // interpolating for visualization
                        psi_node(n) = interpolation_LSQR_Cell_Centered(my_quadtree, psi_n, x, y, NEUMANN, level_set_n);
                        psi_exact_node(n) = interpolation_LSQR_Cell_Centered(my_quadtree, psi_exact, x, y, NEUMANN, level_set_n);

                    }

                char file_name [500];

                sprintf(file_name,  "%s/cell_splitting_%d.vtk",FullPath.c_str(),n);
                my_quadtree.print_VTK_Format(file_name);
                my_quadtree.print_VTK_Format(level_set_n, "level_set",file_name);
                my_quadtree.print_VTK_Format(psi_node, "psi_node",file_name);
                my_quadtree.print_VTK_Format(psi_exact, "psi_exact_node",file_name);

////////////////////////////////////////////////////Error-Analysis////////////////////////////////////////////////////////
                cout<< "+++++++++++++++++++++++++" << endl;
                cout<<"VOLUME: "<< ls_nodes.area_In_Negative_Domain() << endl;
                err_t = ABS(total_psi - total_psi_exact);
//                err_t = ABS(total_psi - total_psi_exact)/pow(2,max_level);

                // this is not working the way I had imagined
//                err_t = ABS(total_psi - total_psi_exact)/ls_nodes.volume_In_Negative_Domain();

                cout<<"Error:"<< err_t <<endl;
                cout<< "+++++++++++++++++++++++++" << endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                }

        filestream << n << "     " << err_t << endl;

        n++;
        t+=dt;
        double smallest_dx = xL/pow(2,max_level);
        dt = dt_const * pow(smallest_dx,2);

        }
    
    return 0;
}

