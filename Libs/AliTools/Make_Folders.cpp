// this will be the main file
// creates folders automatically + textfiles for mass conservations

#include <lib/AliTools/Make_Folders.h>

using namespace std;


void MakeFolder(char hole, char DiffRXN,char DiffOnly,double D_psi, double D_zeta, double gamma_AtoB, double gamma_BtoA, int max_level, int min_level, double tOrder,
                string* txtPathA,string* txtPathB,string* daughterPathA,string* daughterPathB,string* motherPathB,string* motherPathA,char* FolderPath,
                char initCond, string FullPath,string* Return_FolderPAth)
{
    cout<< "Directory Char : h:" << hole << "** D-R:"<< DiffRXN <<
           " ** D: " << DiffOnly << endl;

    if (hole == 'T')
    {
        string Dpsi = "_Dpsi_" + to_string(D_psi);
        string Dzeta = "_Dzeta_" + to_string(D_psi);
        string gammaAB = "_gAB_" + to_string(gamma_AtoB);
        string gamma_BA = "_gBA_" + to_string(gamma_BtoA);
        string max_mesh = "MAX" + to_string(max_level);
        string min_mesh = "_MIN" + to_string(min_level);
        string initial_condition = "_InitCond_";
        string dt = "_dt_" + to_string(tOrder);
        initial_condition += initCond;

        string folder = max_mesh + min_mesh + dt + Dpsi + Dzeta + gammaAB + gamma_BA + initial_condition ;
        string PATH = "/Users/aliheydari/Box/Prion Simulations/Hole/";
        FullPath = PATH + folder;
        *Return_FolderPAth = FullPath;
        //      NAME = FullPath;
        
        FolderPath = const_cast<char*>(FullPath.c_str());
        
        if (mkdir(FolderPath, 0777) == -1)
            cerr << "Error :  " << strerror(errno) << endl;
        
        else
            cout << "***********Directory created***********" << endl
            << "Path: " <<  FolderPath << endl;
        
        // text file for the masses
        
        *txtPathA = PATH + "TotalMassData/" + folder + "_forA.txt";
        *txtPathB = PATH + "TotalMassData/" + folder + "_forB.txt";
        *daughterPathA = PATH + "TotalMassDaughter/" + folder + "_forA.txt";
        *daughterPathB = PATH + "TotalMassDaughter/" + folder + "_forB.txt";
        *motherPathB = PATH + "TotalMassMother/" + folder + "_forA.txt";
        *motherPathA = PATH + "TotalMassMother/" + folder + "_forB.txt";
        
    }
    
    
    else if (DiffRXN == 'T')
    {
        
        string Dpsi = "_Dpsi_" + to_string(D_psi);
        string Dzeta = "_Dzeta_" + to_string(D_zeta);
        string gammaAB = "_gAB_" + to_string(gamma_AtoB);
        string gamma_BA = "_gBA_" + to_string(gamma_BtoA);
        string max_mesh = "MAX" + to_string(max_level);
        string min_mesh = "_MIN" + to_string(min_level);
        string dt = "_dt_" + to_string(tOrder);
        string initial_condition = "_InitCond_";
        initial_condition += initCond;
        
        string folder = max_mesh + min_mesh + dt + Dpsi + Dzeta + gammaAB + gamma_BA + initial_condition ;
        string PATH = "/Users/aliheydari/Box/Prion Simulations/DiffRxn/";
        FullPath = PATH + folder;
        *Return_FolderPAth = FullPath;
        FolderPath = const_cast<char*>(FullPath.c_str());
        
        if (mkdir(FolderPath, 0777) == -1)
            cerr << "Error :  " << strerror(errno) << endl;
        else
            cout << "***********Directory created***********" << endl
            << "Path: " << FolderPath << endl;
        
        
        // text file for the masses
        
        *txtPathA = PATH + "TotalMassData/" + folder + "_forA.txt";
        *txtPathB = PATH + "TotalMassData/" + folder + "_forB.txt";
        *daughterPathA = PATH + "TotalMassDaughter/" + folder + "_forA.txt";
        *daughterPathB = PATH + "TotalMassDaughter/" + folder + "_forB.txt";
        *motherPathB = PATH + "TotalMassMother/" + folder + "_forA.txt";
        *motherPathA = PATH + "TotalMassMother/" + folder + "_forB.txt";
        
        //        return txtPathA,txtPathB,daughterPathA,daughterPathB,motherPathB,motherPathA;
    }
    
    
    else if (DiffOnly == 'T')
    {
        string Dpsi = "_Dpsi_" + to_string(D_psi);
        string Dzeta = "_Dzeta_" + to_string(D_zeta);
        string gammaAB = "_gAB_" + to_string(gamma_AtoB);
        string gamma_BA = "_gBA_" + to_string(gamma_BtoA);
        string max_mesh = "MAX" + to_string(max_level);
        string min_mesh = "_MIN" + to_string(min_level);
        string dt = "_dt_" + to_string(tOrder);
        string initial_condition = "_InitCond_";
        initial_condition += initCond;
        
        string folder = max_mesh + min_mesh + dt + Dpsi + Dzeta + initial_condition ;
        string PATH = "/Users/aliheydari/Box/Prion Simulations/DiffOnly/";
        FullPath = PATH + folder;
        *Return_FolderPAth = FullPath;
        FolderPath = const_cast<char*>(FullPath.c_str());
        
        if (mkdir(FolderPath, 0777) == -1)
            cerr << "Error :  " << strerror(errno) << endl;
        
        else
            cout << "***********Directory created***********" << endl
            << "Path: " << PATH << endl;
        
        
        // text file for the masses
        *txtPathA = PATH + "TotalMassData/" + folder + "_forA.txt";
        *txtPathB = PATH + "TotalMassData/" + folder + "_forB.txt";
        *daughterPathA = PATH + "TotalMassDaughter/" + folder + "_forA.txt";
        *daughterPathB = PATH + "TotalMassDaughter/" + folder + "_forB.txt";
        *motherPathA = PATH + "TotalMassMother/" + folder + "_forA.txt";
        *motherPathB = PATH + "TotalMassMother/" + folder + "_forB.txt";
        
        //        return txtPathA,txtPathB,daughterPathA,daughterPathB,motherPathB,motherPathA;
        
    }
    
    
    else
    {
        cout <<"Not sure which test we are running here..." << endl
        <<" Can NOT create a directory for the output";
        throw invalid_argument("[Make_Folder Error]: Did not provide the system type (Hole,DiffOnly or DiffRxn) correctly \n");
        exit (1);
    }
    
}




