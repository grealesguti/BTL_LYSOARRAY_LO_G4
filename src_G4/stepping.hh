#ifndef STEPPING_HH
#define STEPPING_HH

// Include necessary Geant4 and custom headers
#include "G4UserSteppingAction.hh" // Base class for user stepping actions
#include "G4Step.hh" // Class representing a step in Geant4
#include "G4Track.hh" // Class representing a track in Geant4
#include "G4SystemOfUnits.hh" // System of units for Geant4

// Custom headers for construction, event handling, detector, and particle types
#include "construction.hh"
#include "event.hh"
#include "detector.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4OpticalPhoton.hh"
#include "G4Args.hh" // Custom header for argument handling

// Declare the MySteppingAction class, inheriting from G4UserSteppingAction
// This class is used to define custom actions that occur at each step of a particle's path through the simulation environment.
class MySteppingAction : public G4UserSteppingAction
{
public:
    // Constructor takes pointers to MyEventAction and MyG4Args, for configuration and event handling
    MySteppingAction(MyEventAction* eventAction, MyG4Args*);
    // Destructor
    ~MySteppingAction();

    // Virtual method to perform custom actions at each step of a particle's path
    virtual void UserSteppingAction(const G4Step*);

private:
    // Pointer to MyEventAction for handling event-related actions
    MyEventAction *fEventAction;
    // Pointer to MyG4Args for passing arguments
    MyG4Args* PassArgs;
};

#endif // STEPPING_HH
