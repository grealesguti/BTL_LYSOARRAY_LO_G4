#ifndef PHYSICS_HH
#define PHYSICS_HH

// Include necessary Geant4 headers for physics processes
#include "G4VModularPhysicsList.hh" // Base class for modular physics lists
#include "G4EmStandardPhysics.hh" // Standard electromagnetic physics
#include "G4OpticalPhysics.hh" // Optical physics
#include "G4DecayPhysics.hh" // Decay physics
#include "G4RadioactiveDecayPhysics.hh" // Radioactive decay physics
#include "G4StepLimiterPhysics.hh" // Step limiter physics

// Declare the MyPhysicsList class, inheriting from G4VModularPhysicsList
// This class is used to define a custom physics list for a Geant4 simulation.
class MyPhysicsList : public G4VModularPhysicsList
{
public:
    // Constructor and destructor
    MyPhysicsList();
    ~MyPhysicsList();

    // Methods to add or remove physics processes
    // These methods are not explicitly defined here but are likely inherited from G4VModularPhysicsList
    // and can be overridden to customize the physics processes included in the simulation.
};

#endif // PHYSICS_HH
