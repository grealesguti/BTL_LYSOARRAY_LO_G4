#ifndef ACTION_HH
#define ACTION_HH

// Include necessary Geant4 and custom headers
#include "G4VUserActionInitialization.hh" // Base class for user action initialization
#include "generator.hh" // Custom header for generator actions
#include "run.hh" // Custom header for run actions
#include "event.hh" // Custom header for event actions
#include "stepping.hh" // Custom header for stepping actions
#include "tracking.hh" // Custom header for tracking actions
#include "G4Args.hh" // Custom header for argument handling

// Declare the MyActionInitialization class, inheriting from G4VUserActionInitialization
class MyActionInitialization : public G4VUserActionInitialization
{
public:
    // Constructor takes a pointer to MyG4Args, for configuration (loading input variables)
    MyActionInitialization(MyG4Args*);
    // Destructor
    ~MyActionInitialization();
    
    // Virtual method to build user actions, including primary generator, run, and stepping actions
    virtual void Build() const;

private:
    // String to specify the output name
    G4String OutputName;
    // Pointer to MyG4Args for passing arguments
    MyG4Args* PassArgs;
};

#endif // ACTION_HH
