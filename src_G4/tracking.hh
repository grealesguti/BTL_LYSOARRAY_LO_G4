#ifndef TRACKING_HH
#define TRACKING_HH

// This class is designed to handle tracking actions, 
// which are crucial for monitoring the behavior of 
// particles as they move through the simulation environment.

// Include necessary Geant4 and custom headers
#include "G4UserEventAction.hh" // Base class for user event actions
#include "G4Track.hh" // Class representing a track in Geant4
#include "G4AnalysisManager.hh" // Class for managing analysis tools
#include "G4OpticalPhoton.hh" // Class for optical photons
#include "G4UserTrackingAction.hh" // Base class for user tracking actions
#include "G4TrackingManager.hh" // Class for managing tracking
#include "event.hh" // Custom header for event actions
#include "G4Args.hh" // Custom header for argument handling
#include "run.hh" // Custom header for run actions

// Declare the MyTrackingAction class, inheriting from G4UserTrackingAction
class MyTrackingAction : public G4UserTrackingAction
{
public:
    // Constructor takes pointers to MyEventAction and MyG4Args, for configuration and event handling
    MyTrackingAction(MyEventAction* eventAction, MyG4Args*);
    // Destructor
    ~MyTrackingAction();

    // Virtual method to perform actions before tracking a particle
    virtual void PreUserTrackingAction(const G4Track*);
    // Virtual method to perform actions after tracking a particle
    virtual void PostUserTrackingAction(const G4Track*);

private:
    // Counter for tracking photon counts
    G4double trPhCount;
    // Pointer to MyEventAction for handling event-related actions
    MyEventAction *fEventAction;
    // Pointer to MyG4Args for passing arguments
    MyG4Args* PassArgs;
};

#endif // TRACKING_HH
