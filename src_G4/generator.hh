#ifndef GENERATOR_HH
#define GENERATOR_HH

// Include necessary Geant4 headers for primary particle generation
#include "G4VUserPrimaryGeneratorAction.hh" // Base class for user primary generator actions
#include "G4ParticleGun.hh" // Class for generating primary particles
#include "G4SystemOfUnits.hh" // System of units for Geant4
#include "G4ParticleTable.hh" // Class for particle table
#include "G4Geantino.hh" // Class for Geantino particles
#include "G4IonTable.hh" // Class for ion table
#include "G4GenericMessenger.hh" // Class for handling command-line arguments

// Custom headers for construction and argument handling
#include "construction.hh"
#include "G4Args.hh"

// Declare the MyPrimaryGenerator class, inheriting from G4VUserPrimaryGeneratorAction
// This class is used to define a custom primary particle generator for a Geant4 simulation.
class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    // Constructor takes a pointer to MyG4Args, for configuration
    MyPrimaryGenerator(MyG4Args*);
    // Destructor
    ~MyPrimaryGenerator();
    
    // Virtual method to generate primary particles for a simulation event
    virtual void GeneratePrimaries(G4Event*);

private:
    // Pointer to G4ParticleGun for generating primary particles
    G4ParticleGun *fParticleGun;
    // Position of the particle gun along the z-axis
    G4double zParPos;
    // Pointer to G4GenericMessenger for handling command-line arguments
    G4GenericMessenger *fMessenger;
    // Pointer to MyG4Args for passing arguments
    MyG4Args* PassArgs;
};

#endif // GENERATOR_HH
