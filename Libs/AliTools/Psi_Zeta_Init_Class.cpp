#include<lib/AliTools/Psi_Zeta_Init_Class.h>
//// Space discretization

using namespace CASL;

Initial_Solution::Initial_Solution(){
    lip = 1.;
}


double Initial_Solution::operator() (double x, double y, double z) const
{
    double r1_start = 0.25;

    double x1_start = xL/2.;
    double y1_start = 0.25;
    double z1_start = zL/2.;
    double phi = sqrt(SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))-(r1_start);

    if (phi>0.)
        return 0;
    else
        // Gaussian Dist, exact same as Zeta
        return ((ABS(10* exp(-1 * 100 * (SQR(x-x1_start)+SQR(y-y1_start)+SQR(z-z1_start))))));
}

zeta_Initial_Solution::zeta_Initial_Solution(){
    lip = 1.;
}


double zeta_Initial_Solution::operator() (double x, double y, double z) const
{
    double r1_start = 0.25;

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



