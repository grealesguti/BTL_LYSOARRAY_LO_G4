#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
//#include "construction.hh"
#include "Randomize.hh"
//#include "sim.cc"
#include "G4Args.hh"


#include <string.h>
#include "G4AnalysisManager.hh"
//#include "g4root.hh" // NOT FOUND, G4AnalysisManager include as solution

class MyRunAction : public G4UserRunAction
{
public:
    MyRunAction(G4String ,MyG4Args*);
    ~MyRunAction();
	//void Merge(const G4Run* run) override;

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
private :
    G4String command, OutputName;
    //std::vector<double> eventValues; // data member to store the values for each event

    MyG4Args* PassArgs;
};

#endif
