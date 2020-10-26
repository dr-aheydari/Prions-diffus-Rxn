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
#include <vector>

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
#include <math.h>

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
//double T = 9.5; //final time
//double dt = T/500.; //time step

double T = 9.5; //final time
double dt = T/100.; //time step


double tOrder = dt;
double t = 0.0; //kind of the same as dt
int n = 0; // just for iterations

//// Space discretization
double xL = 0.75; //space [0,L]
double yL = 0.75;
double zL = 0.75;


////////////
// run min 5, max 7 for camera ready

int min_level = 3; //for our adaptive meshing
int max_level = 7; //for our adaptive meshing
////////////

// System parameters
double D_psi = 0.001; // Diffusion coeffcient for "healthy" protein
double D_zeta = 0.001; // Diffusion coeffcient for aggregate

//rate
double gamma_AtoB = 0.001*0 ;
double initial_pop = 10*0;
double gamma_BtoA = 0.01*0;
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


//double direction = 0.1 ; //rand() / (RAND_MAX + 1. + epsilon);
double d_5 = 0.0;
double d_6 = 0.0;

// WHAT IS THIS?
double delta_phi = .01;

/* critical "shell" we draw for our geometry
 if distance between a compartment and the geometry is less than this,
 then we consider it to be like a "collision" and we count it
 */
double epsilon_phi = 0.002;
double k = 1;

double shounter = 0;
double shounter2 = 0;
double hounter = 0;
double hounter2 = 0;
ArrayV <double> velo(13);
int global_collision = 0;
int global_collision2 = 0;

double global_dist_5;
double global_dist_6;



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
    
    //    double b1 = a1;
    double b2 = 1.;
    
    double c1 = 0.9;
    double c2 = 0.9;
    
    double x1_ell = (xL/2.2) + 0.03; //where we center sphere
    double y1_ell = (yL/2.2) - 0.17;
    double z1_ell = (zL/2.2) -0.02 ;
    
    double x2_ell = xL/2.2+0.07;
    double y2_ell = yL/2.2-0.17;
    double z2_ell = zL/2.2+0.07;
    
    
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
    
    double operator() (double x, double y, double z) const
    {
        
        auto [phi1, phi2, phi3,phi4 ,phi5, phi6] =  Assign_phi(x, y,z,velo);
        
        
        if (hole=='F')
        {
            return MIN(phi1,phi2);
        }
        
        else if (hole=='T')
        {
            
            
            return MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi5,-1*phi6);
            
            //            if (compart == 'T')
            //            {
            
            //                return MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi5,-1*phi6);
            
            //            }
            
            //            else if (compart == 'F')
            //            {
            //                return MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4));
            //            }
            
        }
        
        else
        {
            cout<<"Undefined Case of object in the middle...abort" << endl;
            return -1;
        }
        
    }
    
    
    
    tuple<double, double, double, double , double, double > Assign_phi(double x, double y, double z, ArrayV <double> velocity) const
    {
        
        
        double phi1,phi2,phi3,phi4,phi5,phi6;
        
        // (x-x1_start) *(y-y1_start) * (z-z1_start); //mother
        //double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t)); //daughter // sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*t);
        
        // Changed it to y1 start so that we start at the center of the mother cell
        // to make the daughter cell not exit the domain
        
        
        phi1 = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start);
        //        phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start))+SQR(z-z1_start))-(r1_start/3);
        
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
            return {phi1, phi2, phi3,phi4,phi5,phi6};
        }
        
        else if (hole == 'T')
        {
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
            
            // if we have compartments we want to take into account
            if (compart == 'T')
            {
                // count the number of collision for the two compartments
                int counter;
                int bounter;
                // random angles of rotation for the first compartment
                double theta1 = velo(7);
                double omega1 = velo(8);
                double rho1 = velo(9);
                auto [xx_1,xy_1,xz_1,yx_2,yy_2,yz_2,zx_3,zy_3,zz_3] = Rotation_Const(theta1,omega1,rho1);
                
                // random angles of rotation for the second compartment
                double theta2 = velo(10);
                double omega2 = velo(11);
                double rho2 = velo(12);
                
                // these values will store the last spot where we had no issues
                double last_x1;
                double last_y1;
                double last_z1;
                double last_x2;
                double last_y2;
                double last_z2;
                
                if (t < .32)
                {
                    counter = 0;
                    
                    phi5 = sqrt(SQR((x-x1_ell)/a1) +SQR((y-(y1_ell))/b1) + SQR((z-(z1_ell))/c1))-0.015;
                    
                    phi6 = sqrt(SQR((x-(x2_ell))/a1) + SQR((y-(y2_ell))/b1) + SQR((z-(z2_ell))/c1))-0.02;
                    
                    last_x1 = x1_ell;
                    last_y1 = y1_ell;
                    last_z1 = z1_ell;
                    
                    last_x2 = x1_ell;
                    last_y2 = y1_ell;
                    last_z2 = z1_ell;
                    
                    
                }
                
                else
                {
                    // velocity directions for first compartment
                    double x_dir;
                    double y_dir;
                    double z_dir;
                    
                    // velocity compartment for the other compartment
                    double x_dir2;
                    double y_dir2;
                    double z_dir2;
                    
                    int counter = 15;
                    double delta_t = 0.03 ;
                    
                    //                 phi5 = sqrt(SQR((x-x1_ell)/a1) +SQR((y-y1_ell)/b1) + SQR((z-z1_ell)/c1))-0.015;
                    
                    while (counter > 0 && delta_t > 0.001)
                    {
                        
                        x_dir = (1*(velocity(1)*delta_t));
                        y_dir = (1*(velocity(2)*delta_t));
                        z_dir = (1*(velocity(3)*delta_t));
                        
                        
                        //                     phi5 = sqrt(SQR((x-(x1_ell+(x_dir)))/a1)+SQR((y-(y1_ell+(y_dir)))/b1)+
                        //                                                             SQR((z-(z1_ell+(z_dir)))/c1))-0.015;
                        
                        
                        phi5 = sqrt(SQR(((x-(x1_ell+(x_dir)))*xx_1 + (y-(y1_ell+(y_dir)))*xy_1 + (z-(z1_ell+(z_dir)))*xz_1)/a1)
                                    
                                    +SQR(((x-(x1_ell+(x_dir)))*yx_2 + (y-(y1_ell+(y_dir)))*yy_2 + (z-(z1_ell+(z_dir)))*yz_2)/b1)
                                    
                                    +SQR(((x-(x1_ell+(x_dir)))*zx_3 + (y-(y1_ell+(y_dir)))*zy_3 + (z-(z1_ell+(z_dir)))*zz_3)/c1))-0.015;
                        
                        // approximation of distance between the two: add dist + dist2
                        
                        counter = Find_collision(x, y, z ,'5',phi5,phi6);
                        global_dist_5 = Find_distance(x, y, z ,'5',phi5);
                        
                        delta_t = delta_t/2;
                        global_collision = counter;
                        shounter = delta_t;
                        
                    }
                    
                    //                 cout << phi5 << endl;
                    
                    
                    
                    if (counter == 0)
                    {
                        
                        //                     phi5 = sqrt(SQR((x-(x1_ell+(x_dir)))/a1)+SQR((y-(y1_ell+(y_dir)))/b1)+
                        //                                                             SQR((z-(z1_ell+(z_dir)))/c1))-0.015;
                        
                        phi5 = sqrt(SQR(((x-(x1_ell+(x_dir)))*xx_1 + (y-(y1_ell+(y_dir)))*xy_1 + (z-(z1_ell+(z_dir)))*xz_1)/a1)
                                    
                                    +SQR(((x-(x1_ell+(x_dir)))*yx_2 + (y-(y1_ell+(y_dir)))*yy_2 + (z-(z1_ell+(z_dir)))*yz_2)/b1)
                                    
                                    +SQR(((x-(x1_ell+(x_dir)))*zx_3 + (y-(y1_ell+(y_dir)))*zy_3 + (z-(z1_ell+(z_dir)))*zz_3)/c1))-0.015;
                        
                        last_x1 = x1_ell+(1*(velocity(1)*delta_t));
                        last_y1 = y1_ell+(1*(velocity(2)*delta_t));
                        last_z1 = z1_ell+(1*(velocity(3)*delta_t));
                        
                        
                        
                    }
                    else
                    {
                        
                        phi5 = sqrt(SQR((x-last_x1)/a1) + SQR((y-(last_y1))/b1) + SQR((z-last_z1)/c1))-0.015;
                        //                      phi5 = sqrt(SQR((x-x1_ell)/a1) +SQR((y-(y1_ell))/b1) + SQR((z-(z1_ell))/c1))-0.015;
                        
                    }
                    
                    
                    // now for phi_6
                    global_collision = counter;
                    delta_t = 0.04;
                    bounter = 15;
                    
                    // random angles of rotation for the second compartment
                    double theta2 = velo(10);
                    double omega2 = velo(11);
                    double rho2 = velo(12);
                    
                    
                    auto [xx_11,xy_11,xz_11,yx_22,yy_22,yz_22,zx_33,zy_33,zz_33] = Rotation_Const(theta2,omega2,rho2);
                    
                    
                    
                    while (bounter > 0 && delta_t > 0.001)
                    {
                        
                        // for now we use the same velocity vector as before
                        x_dir2 = (1*(velocity(4)*delta_t));
                        y_dir2 = (1*(velocity(5)*delta_t));
                        z_dir2 = (1*(velocity(6)*delta_t));
                        
                        
                        //                     phi6 = sqrt(SQR((x-(x1_ell+0.07))/a1)+SQR((y-(y1_ell-0.17))/b1)+SQR((z-(z1_ell+0.07))/c1))-0.025;
                        
                        //                     phi6 = sqrt(SQR((x-(x2_ell+(x_dir2)))/a1)+SQR((y-(y2_ell+(y_dir2)))/b1)+
                        //                                 SQR((z-(z2_ell+(1*(z_dir2))))/c1))-0.02;
                        
                        
                        phi6 = sqrt(SQR(((x-(x2_ell+(x_dir2)))*xx_11 + (y-(y2_ell+(y_dir2)))*xy_11 + (z-(z2_ell+(z_dir2)))*xz_11)/a1)
                                    
                                    +SQR(((x-(x2_ell+(x_dir2)))*yx_22 + (y-(y2_ell+(y_dir2)))*yy_22 + (z-(z2_ell+(z_dir2)))*yz_22)/b1)
                                    
                                    +SQR(((x-(x2_ell+(x_dir2)))*zx_33 + (y-(y2_ell+(y_dir2)))*zy_33 + (z-(z2_ell+(z_dir2)))*zz_33)/c1))-0.02;
                        
                        
                        
                        
                        // approximation of distance between the two: add dist + dist2
                        
                        //                     counter = Find_collision(x, y, z, counter ,'5',phi5);
                        global_dist_6 = Find_distance(x, y, z ,'6',phi6);
                        bounter = Find_collision(x,y,z,'6',phi5,phi6);
                        
                        delta_t = delta_t/2;
                        shounter2 = delta_t;
                        
                    }
                    
                    if (bounter == 0)
                    {
                        
                        //                     phi6 = sqrt(SQR((x-(last_x2))/a1)+SQR((y-(last_y2))/b1)+SQR((z-(last_z2))/c1))-0.025;
                        
                        //                     phi6 = sqrt(SQR((x-(x2_ell+(x_dir2)))/a1)+SQR((y-(y2_ell+(y_dir2)))/b1)+
                        //                                 SQR((z-(z2_ell+(1*(z_dir2))))/c1))-0.02;
                        
                        phi6 = sqrt(SQR(((x-(x2_ell+(x_dir2)))*xx_11 + (y-(y2_ell+(y_dir2)))*xy_11 + (z-(z2_ell+(z_dir2)))*xz_11)/a1)
                                    
                                    +SQR(((x-(x2_ell+(x_dir2)))*yx_22 + (y-(y2_ell+(y_dir2)))*yy_22 + (z-(z2_ell+(z_dir2)))*yz_22)/b1)
                                    
                                    +SQR(((x-(x2_ell+(x_dir2)))*zx_33 + (y-(y2_ell+(y_dir2)))*zy_33 + (z-(z2_ell+(z_dir2)))*zz_33)/c1))-0.02;
                        
                        
                        
                        last_x2 = x2_ell+(1*(velocity(1)*delta_t));
                        last_y2 = y2_ell+(1*(velocity(2)*delta_t));
                        last_z2 = z2_ell+(1*(velocity(3)*delta_t));
                        
                        
                    }
                    
                    else
                    {
                        phi6 = sqrt(SQR((x-(last_x2))/a1) + SQR((y-(last_y2))/b1) + SQR((z-(last_z2))/c1))-0.02;
                        //                     cout << last_x2;
                    }
                    
                    
                    global_collision2 = bounter;
                    
                }
                
                
                
            }
            
            
            return {phi1, phi2, phi3, phi4 ,phi5, phi6};
        }
        
        
        
        
    }
    
    
    
    
    //    double Assign_desiredPhi(double x, double y, double z, double dir, double forcing_term)
    //    {
    //       double movement = (1+0.15*t) + forcing_term;
    //       double phi5 = sqrt(SQR((x-x1_ell*movement)/a1)+SQR((y-(y1_ell-0.17))/b1)+SQR((z-z1_ell)/c1))-0.02;
    //       return phi5;
    //    }
    
    int Find_Intersect(double x, double y, double z, int counter ,char whichOne,ArrayV <double> direction) const
    {
        auto [phi1, phi2,phi3,phi4,phi5,phi6] = Assign_phi(x,y,z,direction);
        
        if (whichOne == '5')
        {
            // define all the geometry minus the desired level-set function
            
            //
            if (sgn(phi5) == -1 && sgn(MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi6)) == 1)
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
        
        
        if (whichOne == '6')
        {
            // define all the geometry minus the desired level-set function
            
            //
            if (sgn(phi6) == -1 && sgn(MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi5)) == 1)
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
        else
        {
            // the right option was not provided
            cout<<"Wrong input" << endl;
            return 0;
        }
    }
    
    
    
    
    
    tuple<double, double, double, double > Assign_phi4(double x, double y, double z) const
    {
        
        //        double direction;
        
        
        double phi1, phi2,phi3,phi4;
        // (x-x1_start) *(y-y1_start) * (z-z1_start); //mother
        //double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t)); //daughter // sqrt(SQR(x-x2_start)+SQR(y-(y2_start*(1+0.35*t)))+SQR(z-z2_start))-(r2_start+alpha2*t);
        
        // Changed it to y1 start so that we start at the center of the mother cell
        // to make the daughter cell not exit the domain
        
        
        phi1 = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start);
        //        phi3 = sqrt(SQR(x-x1_start)+SQR(y-(y1_start))+SQR(z-z1_start))-(r1_start/3);
        if (t < threshold)
        {
            phi2= sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*t)))+SQR(z-z2_start))-(r2_start+alpha2*(t));
        }
        else
        {
            phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y1_start*(1+0.15*threshold)))+SQR(z-z2_start))-(r2_start+alpha2*(threshold));
            
        }
        //        double phi2 = sqrt(SQR(x-x2_start)+SQR(y-(y2_start))+SQR(z-z2_start))-(r2_start);
        
        if (hole == 'F')
        {
            cout<< " No Holes" << endl;
            return {phi1, phi2,phi3,phi4};
        }
        
        if (hole == 'T')
        {
            
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
            
            return {phi1, phi2,phi3,phi4};
            
        }
        
        
        
        
    }
    
    
    
    
    double Find_distance(double x, double y, double z,char whichOne, double phi_comp) const
    {
        auto [phi1, phi2,phi3,phi4] = Assign_phi4(x,y,z);
        
        // define all the geometry minus the desired level-set function
        double dist1;
        double dist2;
        double dist;
        if (whichOne == '5')
        {
            ArrayV <double> dist_vec(my_octree.number_Of_Nodes());
#pragma omp parallel for
            for(int i = 0; i<my_octree.number_Of_Nodes(); i++)
            {
                dist1 = ABS(phi_comp);
                dist2 = ABS(MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4)));
                dist_vec(i) = ABS(dist1 - dist2);
                
            }
            
            dist = Min_of_Vector(dist_vec, my_octree.number_Of_Nodes());
            int collisions = 0;
            dist = MIN(100.0,dist);
            
            if (dist <= epsilon_phi)
            {
                collisions ++;
                hounter = 0;
            }
            
            else
            {
                // should break and loop out
                hounter = dist;
            }
            
            return dist;
        }
        
        else if (whichOne == '6')
        {
            ArrayV <double> dist_vec(my_octree.number_Of_Nodes());
#pragma omp parallel for
            for(int i = 0; i<my_octree.number_Of_Nodes(); i++)
            {
                double dist1 = ABS(phi_comp);
                double dist2 = ABS(MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4)));
                dist_vec(i) = ABS(dist1 - dist2);
                
            }
            double dist = Min_of_Vector(dist_vec, my_octree.number_Of_Nodes());
            int collisions = 0;
            dist = MIN(100.0,dist);
            
            if (dist <= epsilon_phi)
            {
                collisions ++;
                hounter = 0;
            }
            
            else
            {
                // should break and loop out
                hounter = dist;
            }
            
            return dist;
            
        }
        
    }
    
    
    double Min_of_Vector(ArrayV <double> arr, int n) const
    {
        
        double temp = arr(0);
        
#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            if(temp>arr(i))
            {
                temp=arr(i);
            }
        }
        return temp;
    }
    
    
    
    // to return the appropiate constants for multiplying by x y z
    tuple<double, double, double, double, double, double, double, double, double > Rotation_Const(double theta, double omega, double rho ) const
    {
        double xx_1,xy_1,xz_1,yx_2,yy_2,yz_2,zx_3,zy_3,zz_3;
        
        // first row of the matrix (\times x)
        xx_1 = cos(theta)*cos(omega);
        xy_1 = cos(rho)*sin(theta) - cos(theta)*sin(omega)*sin(rho);
        xz_1 = sin(theta)*sin(rho) + cos(theta)*cos(rho)*sin(omega);
        
        // second row of the matrix (\times y)
        yx_2 = -1* cos(omega)*sin(theta);
        yy_2 = cos(theta)*cos(rho) + sin(theta)*sin(omega)*sin(rho);
        yz_2 = cos(theta)*sin(rho) - cos(rho)*sin(theta)*sin(omega);
        
        // third row of the matix (\times z)
        zx_3 = -1*sin(omega);
        zy_3 = -1 *cos(omega)*sin(rho);
        zz_3 = cos(omega)*cos(rho);
        
        
        return {xx_1,xy_1,xz_1,yx_2,yy_2,yz_2,zx_3,zy_3,zz_3};
        
    }
    
    
    
    int Find_collision(double x, double y, double z,char whichOne, double phi5,double phi6) const
    {
        auto [phi1, phi2,phi3,phi4] = Assign_phi4(x,y,z);
        
        int counter = 0;
        if (whichOne == '5')
        {
            // define all the geometry minus the desired level-set function
            
            
            double phi_total = MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi6);
            //
            if (sgn(phi5) == -1 && sgn(phi_total) == 1)
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
        
        else if (whichOne == '6')
        {
            // define all the geometry minus the desired level-set function
            
            
            double phi_total = MAX(MIN(phi1,phi2),-1*MIN(phi3,phi4),-1*phi5);
            //
            if (sgn(phi6) == -1 && sgn(phi_total) == 1)
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
    MakeFolder(compart, hole,DiffRXN,DiffOnly,D_psi, D_zeta, gamma_AtoB, gamma_BtoA, max_level, min_level, tOrder,
               &txtPathA,&txtPathB,&daughterPathA,&daughterPathB,&motherPathB,&motherPathA, FolderPath,initCond,FullPath,&Return_FolderPAth);
    
    FullPath = Return_FolderPAth;
    
    
    char whichOne1= '5';
    char whichOne2= '6';
    
    int coll_counter1 , coll_counter2;
    
    //    get the velocity vector;
    
    
    
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
        daughter_cell(i) = level_set_daughter(x,y,z);
    }
    
    // for loop over all the nodes
    
    ArrayV <double> x_grad(my_octree.number_Of_Nodes());
    ArrayV <double> y_grad(my_octree.number_Of_Nodes());
    ArrayV <double> z_grad(my_octree.number_Of_Nodes());
    
    
    ///////////////////////////////////////////////////////////////////
    //    for (int i = 0; i<my_octree.number_Of_Nodes(); i++)
    //    {
    //       // neighborhood of a node
    //       OctNgbdNodesOfNode nghb;
    //       my_octree.get_Ngbd_Nodes_Of_Node(i,nghb);
    //       double x = my_octree.x_fr_i(my_octree.get_Node(i).i);
    //       double y = my_octree.y_fr_j(my_octree.get_Node(i).j);
    //       double z = my_octree.z_fr_k(my_octree.get_Node(i).k);
    
    //       x_grad(i) = nghb.dx_Central(level_set_n);
    //       y_grad(i) = nghb.dy_Central(level_set_n);
    //       z_grad(i) = nghb.dy_Central(level_set_n);
    //       double forcing_term = k/level_set_n(i) + x_grad(i)+y_grad(i)+z_grad(i);
    ////       cout<< level_set.Assign_desiredPhi(x, y, z, 1, forcing_term) << endl;;
    //    }
    
    //    cout<< z_grad << endl;
    ///////////////////////////////////////////////////////////////////
    
    
    
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
    
    
    
    
    velo.resize_Without_Copy(13);
    //    std::cout << sin(3.14/2) << std::endl;
    // MAIN WHILE LOOP
    while (t <= T)
    {
        
        // generate a velocity vector for the first compartment
        velo(1) = rand() / (RAND_MAX + 1. + epsilon);
        velo(2) = rand() / (RAND_MAX + 1. + epsilon);
        velo(3) = rand() / (RAND_MAX + 1. + epsilon);
        // generate a velocity vector for the second compartment
        velo(4) = rand() / (RAND_MAX + 1. + epsilon);
        velo(5) = rand() / (RAND_MAX + 1. + epsilon);
        velo(6) = rand() / (RAND_MAX + 1. + epsilon);
        
        // generate 3 random roation angles for the first comp.
        velo(7) = rand()  % 101;
        velo(8) = rand()  % 101;
        velo(9) = rand()  % 101;
        
        // generate 3 random roation angles for the second comp.
        velo(10) = rand()  % 101;
        velo(11) = rand()  % 101;
        velo(12) = rand()  % 101;
        
        if (n > 0)
        {
            cout << " \n \n" << endl;
        }
        
        cout<<"We are at step: "<< n <<" at time = "<< t <<endl;
        cout << "-------------Phi5--------------" << endl;
        cout << "Angles of Rotation: (" << velo(7) <<", " << velo(8) <<", "<<velo(10)<<")"<< endl;
        cout << "direction 1 is: (" << velo(1) <<", " << velo(2) <<", "<<velo(3)<<")"<< endl;
        cout << "last min distance was "<< global_dist_5 <<endl;
        cout << "last global collision for phi5 counter: " << global_collision << endl;
        cout << "last dt = " << shounter << endl;
        cout << "+++++++++++++Phi6++++++++++++++" << endl;
        cout << "Angles of Rotation: (" << velo(10) <<", " << velo(11) <<", "<<velo(12)<<")"<< endl;
        cout << "direction 2 is: (" << velo(4) <<", " << velo(5) <<", "<<velo(6)<<")"<< endl;
        cout << "last min distance was "<< global_dist_6 <<endl;
        cout << "last global collision for phi6 counter: " << global_collision2 << endl;
        cout << "last dt = " << shounter2 << endl;
        cout << "=======================================" << endl;
        
        
        
        OcTreeSolverCellCenteredPoisson poisson_solver;
        poisson_solver.set_Octree(my_octree);
        BoundaryConditions3D bc_psi;
        bc_psi.setWallTypes(wall_psi_neumann_type);
        bc_psi.setWallValues(wall_psi_neumann_value);
        bc_psi.setInterfaceType(NEUMANN);
        bc_psi.setInterfaceValue(int_psi_neumann_value);
        
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
            
            
            /////////////////////////////////
            
            // HERE IS WHERE WE PRINT THE NAN checks
            
            //            rhs.CHK_NAN();
            //            cout<<rhs.max_Abs()<<endl;
            //            cout<<rhs.min()<<endl;
            
            /////////////////////////////////
            
            
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
        
        
        // We don't need to output any of these for now!!
        
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
            
            //            // check if the compartments are touching as we update the level-set
            //            coll_counter1 = level_set.Find_Intersect(x,y,z,coll_counter1,whichOne1,velo);
            //            coll_counter2 = level_set.Find_Intersect(x,y,z,coll_counter2,whichOne2,velo);
            
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
        
        ArrayV<double> levels(my_octree.number_Of_Leaves());
#ifdef CASL_OPENMP
#pragma omp parallel for
#endif
        for(CaslInt l=0; l<my_octree.number_Of_Leaves(); l++)
        {
            levels(l) = log( (double) OcTreeStatic::MAX_NUMBER_OF_NODES_IN_ONE_DIRECTION_OCTREE /my_octree.get_Leaf(l).size())/log(2.);
        }
        
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
        my_octree.print_VTK_Format_Cell(levels,"level",file_name);
        
        n++;
        
    }
    
    return 0;
}




