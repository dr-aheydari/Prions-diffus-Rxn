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

#include <lib/arrays/ArrayV.h>
#include <lib/amr/OcTree.h>
#include <lib/amr/OcTreeCellBased.h>
#include <lib/amr/OcTreeLevelSet.h>
#include <lib/amr/OcTreeAdvection.h>
#include <lib/tools/OcTreeToolbox.h>
#include <lib/amr/OcTreeCellBasedLevelSet.h>
#include <lib/solvers/OcTreeSolverCellCenteredPoisson.h>

using namespace std;
using namespace CASL;

// prototype functions :
double Mass_Correction_Diff(double initial_mass,double current_mass);
double SoftDist_Mass_Correction(double initial_mass,double current_mass,double area);



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
double T = 10; //final time
double dt = T/100.; //time step
double t = 0.; //kind of the same as dt
int n = 0; // just for iterations

// Space discretization
double xL = 0.75; //space [0,L]
double yL = 0.75;
double zL = 0.75;
int min_level = 1; //for our adaptive meshing
int max_level = 6; //for our adaptive meshing

// System parameters
double D_psi = 0.001; // Diffusion coeffcient for "healthy" protein
double D_zeta = 0.001; // Diffusion coeffcient for aggregate

//rate
double gamma_var = 0.5;
double initial_pop = 10*0;
double conversion_rate1 = 0.05;
double mu = 0.2 * 0;
//double D_z = 0.001;

double init_mass_A;
double init_mass_B;



// Level_set

class LS : public CF_3
{
public:
    LS() //this is the Level Set Function.
    {
        lip = 1.; //Lipschitz constant, not important for us
    }
    double operator() (double x, double y, double z) const
    {

        double beta = 1.7; //radius of mother and daughter depend on this


        double r1_start = 0.15; //radius

        // so that we have mass conservation:
        double r2_start = 0.05; // + (0.005 * t); //another radius

        double alpha1 = r1_start/T*(beta-1); //placeholders
        double alpha2 = r2_start/T*(beta-1);

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
        if (t < 3)
        {
             phi2= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t));
        }
        else
        {
             phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.35*3)))+SQR(z-z2_start))-(r2_start+alpha2*(t));

        }
        //        double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start))+SQR(z-z2_start))-(r2_start);

        return MIN(phi1,phi2); //this is where the level set comes in, where we're actually doing the level set: we don't have negatives, so whichever gives us zero is our level set (because o(x) = 0, from picture)



    }
} level_set;

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

// Initial solution

class Initial_Solution : public CF_3
{
public:
    Initial_Solution() //initial condition, is a distribution with some concentration (number of species)
    {
        lip = 1.;
    }
    double operator() (double x, double y, double z) const
    {
        double r1_start = 0.15;

        double x1_start = xL/2.;
        double y1_start = 0.25;
        double z1_start = zL/2.;
        double phi = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start);

        if (phi>0.)
            return 0;
        else
            //return ((ABS(1000 * sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))))); //this is where we make it a circle

            // Gaussian Dist, exact same as Zeta
            return ((ABS(10* exp(-1 * 100 * (SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))))));
    }
} psi_initial;


// ZETA INITIAL CONDITION
class zeta_Initial_Solution : public CF_3
{
public:
    zeta_Initial_Solution() //initial condition, is a distribution with some concentration (number of species)
    {
        lip = 1.;
    }
    double operator() (double x, double y, double z) const
    {
        double r1_start = 0.15;

        double x1_start = xL/2.;
        double y1_start = 0.25;
        double z1_start = zL/2.;
        double phi = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start);
        if (phi>0.)
            return 0;
        else
            //return ((ABS(1000 * sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start)))); //this is where we make it a circle

            return ((ABS(10 * exp(-1 * 100 * (SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start)))))); // Gaussian Distribution


    }
} zeta_initial;
















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


// Structs for weight ditribution
struct SoftDist {
    double zeta_add, psi_add,total_mass;
};

//typedef struct greaterSmaller SoftDist;

SoftDist softdist(double initial_totalMass ,double current_A, double current_B,double num_funcs)
{

    // for now we are not gonna mess with num_funcs since we assume its always 2

    SoftDist s;
    double beta = 2;
    double denom = exp(beta* current_A) + exp(beta *current_B);

    double coeff1 = exp(current_A * beta) / (denom) ;
    double coeff2 = exp(current_B * beta) / (denom);

    s.psi_add = coeff1 * initial_totalMass;
    s.zeta_add = coeff2 * initial_totalMass;

    s.total_mass =  coeff1 * initial_totalMass + coeff2 * initial_totalMass;

    return s;
}
















//everything we had before was just declaring classes, but we now need to initialize them

int main()
{
#ifdef CASL_OPENMP
    (void) omp_set_dynamic(false);
    //    if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
    (void) omp_set_num_threads(4);
    cout << "OpenMP - Max number of threads : " << omp_get_max_threads() << endl;
#endif





    my_octree.set_Grid(0.,xL,0.,yL,0.,zL); //setting the grid
    my_octree.construct_Octree_From_Level_Function(level_set, min_level, max_level);
    //    my_octree.construct_Uniform_Octree(max_level);
    my_octree.initialize_Neighbors(); //talking about numerical cells as neighbors of each other


    ArrayV <double>  level_set_n (my_octree.number_Of_Nodes());
    psi_n.resize_Without_Copy(my_octree.number_Of_Leaves());
    psi_nm1.resize_Without_Copy(my_octree.number_Of_Leaves());


    //ZETA
    zeta_n.resize_Without_Copy(my_octree.number_Of_Leaves());
    zeta_nm1.resize_Without_Copy(my_octree.number_Of_Leaves());

#pragma omp parallel for
    for (CaslInt i = 0; i<my_octree.number_Of_Leaves(); i++) //where adaptive meshing is done
    {
        double x = my_octree.x_fr_i(my_octree.get_Leaf(i).icenter());
        double y = my_octree.y_fr_j(my_octree.get_Leaf(i).jcenter());
        double z = my_octree.z_fr_k(my_octree.get_Leaf(i).kcenter());

        // used to be psi_initial, but kist to ,ake srue they have the same initial condition!!

        psi_n(i) = zeta_initial(x,y,z); //first element of psi_n is the initial condition
        //psi_nm1(i) = psi_n(i); //i element of psi_nm1 is i element of psi_n (psi_n is the n+1)

        //ZETA

        // this used to be zeta_initial... but just to make sure it's the same as psi
        zeta_n(i) = zeta_initial(x,y,z);
        //zeta_nm1(i) = 2*gamma_var*zeta_n(i);

    }
    char file_name [500];



    // THIS IS WHAT I ADDED TO DO INTEGRATION OVER THE DOMAIN

    OcTreeLevelSet lev_set;
    lev_set.set_Octree(my_octree);
    //lev_set.set_Phi(psi_n);






    level_set_n.resize_Without_Copy(my_octree.number_Of_Nodes());
#pragma omp parallel for
    for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
    {
        double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
        double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
        double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
        level_set_n(i) = level_set(x,y,z);
    }


    // print to the file for figures
        ofstream filestream("/Users/aliheydari/Documents/Prions Project Data/total_psi_data_diffOnly/A_SoftDist_DStays2.txt");
        ofstream filestream2("/Users/aliheydari/Documents/Prions Project Data/total_zeta_data_diffOnly/B_SoftDist_DStays2.txt");
        ofstream OGstream1("/Users/aliheydari/Documents/Prions Project Data/total_psi_data_diffOnly/Original_SoftDist_DStays2.txt");
        ofstream OGstream2("/Users/aliheydari/Documents/Prions Project Data/total_zeta_data_diffOnly/Original_SoftDist_DStays2.txt");



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



        ArrayV<double> psi_updated;
        psi_updated.resize_Without_Copy(my_octree.number_Of_Leaves());

        //#pragma omp parallel for



                for (int i = 0; i < psi_updated.size(); i++) //because we got adding array operator, just inbetween step
                {
                    // rhs
                    // for now we are doing Diffusion Only
                    psi_updated(i) =psi_n(i) + (initial_pop*dt) - (dt*conversion_rate1*(psi_n(i))*(psi_n(i))) + (2 *dt* gamma_var * zeta_n(i));// MOST RECENT -> DIFFUSION -> psi_n(i) + 2*gamma_var * zeta_nm1(i); // (alpha - beta*(psi_n(i))*(psi_n(i)) + gamma_var * zeta_m(i));//2 * gamma_var * zeta_m(i))


                }


        ArrayV<double> rhs = psi_updated;//previous right hand side, kind of like a placeholder
         //ArrayV<double> rhs = psi_nm1;
        poisson_solver.set_bc(bc_psi);
        poisson_solver.set_Rhs(rhs);
        // the diagonals here are what need to be changed if we want to add reaction terms
        poisson_solver.set_Diagonal_Increment(1.0 + mu * dt);

//        poisson_solver.set_Diagonal_Increment(1.0);

        poisson_solver.set_mu(dt*D_psi);

        poisson_solver.set_Phi(level_set_n);
        poisson_solver.solve(psi_n); //we then solve for psi_n




        OcTreeCellBasedLevelSet ls;
        ls.set_Octree(my_octree);
        ls.set_Phi(level_set_n);
        ls.extrapolate_Along_Normal_Using_Cells(psi_n, bc_psi);
        // defining new variable for the mass conservation
        double total_psi=ls.integral_Cell_Based(psi_n);

        OcTreeLevelSet lvl;
        lvl.set_Octree(my_octree);
        lvl.set_Phi(level_set_n);
        lvl.volume_In_Negative_Domain();



        //ZETA
//        BoundaryConditions3D bc_zeta;
//        bc_zeta.setWallTypes(wall_zeta_neumann_type);
//        bc_zeta.setWallValues(wall_zeta_neumann_value);
//        bc_zeta.setInterfaceType(NEUMANN);
//        bc_zeta.setInterfaceValue(int_zeta_neumann_value);


        // to MAKE SURE EVERYTHING IS EXACTLY THE SAME

        BoundaryConditions3D bc_zeta;
        bc_zeta.setWallTypes(wall_psi_neumann_type);
        bc_zeta.setWallValues(wall_psi_neumann_value);
        bc_zeta.setInterfaceType(NEUMANN);
        bc_zeta.setInterfaceValue(int_psi_neumann_value);

        ArrayV<double> zeta_updated;
        zeta_updated.resize_Without_Copy(my_octree.number_Of_Leaves());

        //#pragma omp parallel for



                for (int i = 0; i < psi_updated.size(); i++) //because we got adding array operator, just inbetween step
                {
                    // rhs
                    // since zeta is a dimer, then two monomers give you one dimer so that's why there's the 0.5 in there
                    zeta_updated(i) = zeta_n(i) + 0.5 * (dt*conversion_rate1*(psi_n(i))*(psi_n(i))); // DIFUSION with gamma=0-> zeta_n(i) - 2*gamma_var * zeta_nm1(i); // (alpha - beta*(psi_n(i))*(psi_n(i)) + gamma_var * zeta_m(i));//2 * gamma_var * zeta_m(i))


                }



        //ArrayV<double> rhs_z = zeta_nm1;

        ArrayV<double> rhs_z =zeta_updated;
        poisson_solver.set_bc(bc_zeta);
        poisson_solver.set_Rhs(rhs_z);
        // the diagonals here are what need to be changed if we want to add reaction terms
        poisson_solver.set_Diagonal_Increment(1.0 + (mu * dt) + (dt * gamma_var));

//        poisson_solver.set_Diagonal_Increment(1.);

        poisson_solver.set_mu(dt*D_zeta);

        poisson_solver.set_Phi(level_set_n);
        poisson_solver.solve(zeta_n);

        ls.extrapolate_Along_Normal_Using_Cells(zeta_n, bc_zeta);

        double total_zeta = ls.integral_Cell_Based(zeta_n);

        cout<< "Total Psi (cell based) -> "<< total_psi << endl;
        cout<< "Total zeta (cell based) -> "<< total_zeta << endl;



        if (n == 0)
        {
           total_zeta = total_psi;
           init_mass_A = total_psi;
           cout <<"Initial 'A' Mass set to : " << init_mass_A << endl;
        }


        if (n == 0)
        {
           init_mass_B = total_psi;
           cout <<"Initial 'B' Mass set to : " << init_mass_B << endl;
        }


        //double volume = lvl.volume_In_Negative_Domain();

        //cout << "Volume is " << volume << endl;

        //double corrected_mass = Mass_Correction_Diff(init_mass_A,total_psi,volume);

//        SoftDist dist;

//        dist = softdist(init_mass_A + init_mass_B, total_psi,total_zeta,2);


//        cout << "corrected mass for A is : " << dist.psi_add << endl;


//         //filestream << n << "     " << total_psi << endl;




//        // file output

//         filestream << n << "     " << dist.psi_add << endl;

//         cout << "corrected mass for B is : " << dist.zeta_add << endl;

//         //filestream2 << n << "     " << total_zeta << endl;

//         filestream2 << n << "     " << dist.psi_add << endl;

//        cout << "Total mass is " << dist.total_mass << endl;
//        OGstream1 << n << "     " << total_psi << endl;
//        OGstream2 << n << "     " << total_zeta << endl;







        //SoftDist dist;
        double psi_add = Mass_Correction_Diff(init_mass_A,total_psi);
        //dist = softdist(init_mass_A + init_mass_B, total_psi,total_zeta,2);


        cout << "corrected mass for A is : " << psi_add << endl;

        filestream << n << "     " << psi_add << endl;

        double zeta_add = Mass_Correction_Diff(init_mass_B,total_zeta);
        cout << "corrected mass for B is : " << zeta_add << endl;

        //filestream2 << n << "     " << total_zeta << endl;

        filestream2 << n << "     " << psi_add << endl;

        double total_mass = psi_add + zeta_add;

        cout << "Total mass is " << total_mass << endl;
        OGstream1 << n << "     " << total_psi << endl;
        OGstream2 << n << "     " << total_zeta << endl;





















        t += dt;




        new_octree.make_It_Have_Root_Cell_Only();
        new_octree.set_Grid(0.,xL,0.,yL,0.,zL);
        new_octree.construct_Octree_From_Level_Function(level_set, min_level, max_level);
        //        SplitCriteria criterion;
        //        new_octree.construct_Octree_From_Threshold_Function(criterion, min_level, max_level);


        //        //        new_octree.impose_Uniform_Grid_On_Interface(level_set,max_level,3);
        new_octree.initialize_Neighbors();


        //        SplitCriteria criterion;
        //        new_octree = my_octree;
        //        new_octree.update_Octree(criterion);
        //        new_octree.initialize_Neighbors();


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

            new_psi_n(leaf) = interpolation_LSQR_Cell_Centered(my_octree, psi_n, x, y, z);
            //ZETA
            new_zeta_n(leaf) = interpolation_LSQR_Cell_Centered(my_octree, zeta_n, x, y, z);
        }
        // for visu only
        my_octree = new_octree;
        psi_n.resize_Without_Copy(my_octree.number_Of_Leaves());
        psi_nm1 = new_psi_n;
        psi_n = new_psi_n;
        //ZETA
        zeta_n.resize_Without_Copy(my_octree.number_Of_Leaves());
        zeta_nm1 = new_zeta_n;
        zeta_n = new_zeta_n;



        ArrayV<double> psi_node(my_octree.number_Of_Nodes());
        //ZETA
        ArrayV<double> zeta_node(my_octree.number_Of_Nodes());
#pragma omp parallel for
        for (CaslInt n = 0; n<new_octree.number_Of_Nodes(); n++)
        {
            double x = new_octree.x_fr_i(new_octree.get_Node(n).i);
            double y = new_octree.y_fr_j(new_octree.get_Node(n).j);
            double z = new_octree.z_fr_k(new_octree.get_Node(n).k);

            psi_node(n) = interpolation_LSQR_Cell_Centered(my_octree, psi_n, x, y, z);
            //ZETA
            zeta_node(n) = interpolation_LSQR_Cell_Centered(my_octree, zeta_n, x, y, x);
        }

        level_set_n.resize_Without_Copy(my_octree.number_Of_Nodes());
#pragma omp parallel for
        for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
        {
            double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
            double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
            double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
            level_set_n(i) = level_set(x,y,z);
        }


        // print everything
        sprintf(file_name, "/Users/aliheydari/Documents/Prions Project Data/SoftDist_DStays/cell_splitting_%d.vtk", n);
        my_octree.print_VTK_Format(file_name);
        my_octree.print_VTK_Format(level_set_n, "level_set",file_name);
        my_octree.print_VTK_Format(psi_node, "psi_node",file_name);
        //ZETA
        my_octree.print_VTK_Format(zeta_node, "zeta_node",file_name);
        my_octree.print_VTK_Format_Switch_To_Cell(file_name);
        my_octree.print_VTK_Format_Cell(psi_n, "psi",file_name);
        //ZETA
        my_octree.print_VTK_Format_Cell(zeta_n, "zeta",file_name);


        // print everything

        n++;
    }










    return 0;
}


double Mass_Correction_Diff(double initial_mass,double current_mass)
{

    double epsilon = (initial_mass - current_mass);
    cout << "EPSILON : " << epsilon << endl;
    double updated_mass = (current_mass / (current_mass + epsilon)) * initial_mass;


    return  updated_mass;


}


//double SoftDist_Mass_Correction(double initial_totalMass ,double current_A, double current_B,double num_funcs)
//{
//    double denom = exp(current_A) + exp(current_B);





//    return -1;

//}














