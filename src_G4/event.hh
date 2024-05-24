#ifndef EVENT_HH
#define EVENT_HH

// Include necessary Geant4 headers for event actions and analysis
#include "G4UserEventAction.hh" // Base class for user event actions
#include "G4Event.hh" // Class representing an event in Geant4
#include "G4AnalysisManager.hh" // Class for managing analysis of simulation data
#include "G4UImanager.hh" // Class for user interface management
#include "G4SystemOfUnits.hh" // System of units for Geant4
#include "G4PhysicsOrderedFreeVector.hh" // Class for ordered free vector of physics processes

// Include custom headers for run and construction
#include "run.hh" // Custom header for run-related actions
#include "construction.hh" // Custom header for construction of the simulation world
#include <string.h> // Functions for manipulating C-style strings
#include "G4Args.hh" // Custom header for argument handling

// Declare the MyEventAction class, inheriting from G4UserEventAction
// This class is used to define custom actions to be taken at the beginning and end of each event in a Geant4 simulation.
class MyEventAction : public G4UserEventAction
{
public:
    // Constructor takes pointers to MyRunAction and MyG4Args, presumably for configuration and passing arguments
    MyEventAction(MyRunAction*, MyG4Args*);
    // Destructor
    ~MyEventAction();

    // Virtual methods to define actions at the beginning and end of each event
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    // Methods to add energy deposit, photon count, and linear optical (LO) values
    void AddEdep(G4double edep){fEdep += edep;}
    void AddPh(G4double PhCount){fPhCount += PhCount;}
    void AddLO(G4double LOc){fLO += LOc;}

private:
    // Variables to store energy deposit, photon count, and linear optical (LO) values
    G4double fEdep; // accumulation of energy deposition
    G4double fPhCount; // accumulation of photon counts
    G4double fLO, GenX=0., GenZ=0.; // Additional variables for LO and generic X and Z positions, initialized to 0
    G4double PDE420; // Variable for storing a specific value, possibly related to physics processes (peak PDE value to scale the PDE values from the txt)
    G4PhysicsOrderedFreeVector *PDE; // Pointer to a vector of physics processes for the PDE of the SiPM (read it out)
    G4double PXd, PZd; // Variables for storing generic X and Z positions
    G4String PX, PZ; // Strings for storing generic X and Z positions 
    G4String GLUE_Lstr; // String for storing a glue layer string, possibly related to construction or analysis
    G4String command; // String for storing a command, possibly for analysis or configuration
    MyG4Args* PassArgs; // Pointer to MyG4Args for passing arguments
};

#endif
