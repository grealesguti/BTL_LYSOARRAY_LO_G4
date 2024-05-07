#ifndef G4SIM_HH
#define G4SIM_HH

// Include necessary standard and Geant4 headers
#include <iostream> // For input/output operations
#include "G4RunManager.hh" // Class for managing the simulation run
#include "G4UImanager.hh" // Class for user interface management

// Include custom headers for construction, physics, and action
#include "construction.hh"
#include "physics.hh"
#include "action.hh"
#include <cstdlib> // For standard library functions
#include "G4Args.hh" // Custom header for argument handling
#include <unistd.h> // For exit() function

// Declare the G4simulationNOVIS class
// This class defines each run without loading the visualization modules. Prepared to run in HPCs.
class G4simulationNOVIS 
{
public:
    // Constructor takes command-line arguments, and optional parameters for node configuration
    G4simulationNOVIS(int, char**, G4int Onode = 5, G4int Znode = 1, G4double* radp = NULL);
    // Destructor
    ~G4simulationNOVIS();

    // Method to get the average linear optical (LO) value
    G4double GetAvgLO(){}

    // Methods to get average and standard deviation of linear optical (LO) and linear density (LD) values
    double GetLO_avg(int runid){return ArgInp->GetLOAvg(runid);}
    double GetLO_std(int runid){return ArgInp->GetLOStd(runid);}
    double GetLD_avg(int runid){return ArgInp->GetLDAvg(runid);}
    double GetLD_std(int runid){return ArgInp->GetLDStd(runid);}

private:
    // Pointer to MyG4Args for passing arguments
    MyG4Args *ArgInp;
    // Average linear density (LD) value
    G4double LDavg;
};    

#endif
