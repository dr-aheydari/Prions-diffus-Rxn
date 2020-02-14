////// this will be the header file
////// NaN finder for the values of the total masses
#ifndef PSI_ZETA_INIT_CLASS
#define PSI_ZETA_INIT_CLASS

#include<lib/tools/CASL_math.h>

using namespace CASL;


class Initial_Solution : public CF_3
{
public:
    double xL = 0.75; //space [0,L]
    double yL = 0.75;
    double zL = 0.75;
    Initial_Solution(); //initial condition, is a distribution with some concentration (number of species)
    double operator() (double x, double y, double z) const;
};


// ZETA INITIAL CONDITION
class zeta_Initial_Solution : public CF_3
{
public:
    double xL = 0.75; //space [0,L]
    double yL = 0.75;
    double zL = 0.75;
    zeta_Initial_Solution(); //initial condition, is a distribution with some concentration (number of species)
    double operator() (double x, double y, double z) const ;
};

#endif // PSI_ZETA_INIT_CLASS
