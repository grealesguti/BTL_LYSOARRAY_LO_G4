#ifndef G4ARGS_HH
#define G4ARGS_HH

// Include necessary Geant4 headers for managing the simulation run
#include "G4RunManager.hh" // Class for managing the simulation run

// Include standard C++ headers
#include <cstdlib> // Standard library functions, including dynamic memory management, random number generation, and others
#include <iostream> // Standard input/output stream objects
#include <algorithm> // Template functions for algorithms like std::sort
#include <vector> // Template class for dynamic arrays
#include <fstream> // File stream operations
#include <cstring> // Functions for manipulating C-style strings
#include <ctime> // Functions for time and date
#include <string.h> // Functions for manipulating strings
#include <sstream> // String stream operations
#include <stdio.h> // Standard input/output functions

// Include custom utility header
#include "util.hh" // Custom utility functions and definitions

// Include Geant4 header for user interface management
#include "G4UImanager.hh" // Class for managing the user interface

// Include Geant4 header for logical volumes
#include "G4LogicalVolume.hh" // Class representing a logical volume in Geant4


// This class provides all tuning or modification potential within all other classes.
// This class contains functions to store values within each part of each run and or event of 
// GEANT4. Furthermore, it reads and stores all values introduced by cmd.
// Possible inputs and default values or possible options:
/*
-GeomConfig <int> : options [1 2 3 11 13] This argument changes the configuration of the detector. Configuration 1-3 provide use the default geometry construction of G4 while 11 and 13 use Gmsh and allow shape changes through splines and the -YPos argument. Config. 3 and 13 use 16 crystals while the rest use a single one.
-Zelem <int> : number of mesh nodes along half of Z of the LYSO crystal
-Znode <int> : number of sections of the LYSO mesh, number of <int> within the -Ypos argument minus 1
-Ypos {-<int>-<int>} : value that multiplies the 3mm, by default height of the LYSO, of each spline control node of the mesh.
-incrV <int> : range [1-199] value between 1 and 199 that provides the 1D study of 2 half tetrahedra for the crystal shape decreasing and increasing the LYSO thickness in the middle and edges of the crystal

-o <string> :  name identifies the name.root file created in the Results folder with the information of the simulation.
-runevt <int> : argument that provide how many events will be run. If not given and if you use G4vis.cc you can use the graphical interface.

-rnd <int> <int> <int> : taking 3 arguments in the form of 1 or 0. The second argument sets no random position to the particle impacts, the second sets no random geometry tolerances during the detector construction.
-rndangle : the particle impacts in random directions and not only along the Z direction
-gunmesh <int> <int> : provides a uniform impact pattern in a quarter of the crystal along X and Z. With the given defaults it would run 60 events.
-Muon : sets the particle impact to 2MeV Muons.

-ESRbackpainted : changes the surface coating model so that there is an air gap between the ESR and the LYSO
-ESRdefaultmodel: changes the default G4 model with a polished LYSO, air gap and ESR (LUT DAVIS)
-noESR : no ESR coating is applied to the crystal

-LYSO_L <int> : default 28.5 where the number given is half the lenght of the LYSO crystal along Z
-Acte : command that maintains the XxZ area constant increasint X if LYSO_L is changed
-matchSiPM : command that sets the SiPM width X equal to the LYSO crystal width
-incrSiPM  <int> :range [25-200], increases in percentage the height of the SiPM from its default 3mm 
-SiPM_Z <int>
-RESIN_Z <int>
-RESIN_Z1000 <int>
-Glue_Z <int>
-Glue_Z1000 <int>
-RESIN_W <int>

-GmshView: opens Gmsh before running G4. Right now does not allow to continue the simulation.
-SaveSTL: saves an .stl file of the mesh for further postprocessing

-BC408: changes the crystal material from LYSO to BC408 (polymer).

-TileV0 : changes configuration from bar to tile. (TODO: right now the design is assymmetrical, the triangulation needs to be changed in Gmsh)
-ForceBottomLine: To always be used with -TileV0. Forces a given flat bottom surface.
-SZloc <double> : range [0-1] changes the SiPM location along Z at the bottom of the tile as a percentage of LYSO_L.
*/
class MyG4Args 
{

public:
    MyG4Args(int, char**);
    ~MyG4Args();

    void SetNSGAII();
    
    //*** Initialization functions ***// 
    // -> Set to zero at the begining of run.cc or event.cc to avoid numerical issues or bad counts
    // -> or, define default dimensions
    void InitAllCount(){ArgLO = 0;ArgCrossTalk = 0;TotPh = 0;PhHit=0;Edep=0.;MuonEdep=0.;nPhotL=0;nPhotR=0;PhotTiming[1]=0.;MuonLYSOTrackLength=0.;} // Initialize counters (Energy deposition, photon impacts ... ) to zero for each run
    void InitTotPh(){TotPh = 0;}
    void InitLO(){ArgLO = 0;}
    void InitCT(){ArgCrossTalk = 0;}
    void InitVolume(){Volume = 0;} // Init. counter for LYSO volume (sum of tetrahedron volumes)
    void InitfScoringVolumeVec(std::vector<G4LogicalVolume*> fScoringVolumeVecinit) {fScoringVolumeVec=fScoringVolumeVecinit;}  // Init. scoring volume reference
    void SetLYSOVolumeXY();
    void SetCoordVect();
    void SetYVect(G4double*);
    
    //*** Getters: Used in follow up classes ***//
    // String Name Getters -> used for string creation for output files in following classes ().
    G4String GetOutName() const {return OutName;} 	// Get output name
    G4String GetMacName() const {return MacName;}	// Get macro name
    G4double GetVov() const {return Detection[0];}	// Get overvoltage value
    G4String GetRootFolder() const {return rootfolder;} // Get name of folder to store root files
    
    // Int construction property Getters : true or false properties or integer values
    G4int GetGeomConfig() const {return GeomConfig;}	// Get Geometry configuration. 1 - Single LYSO fully G4 created; 11 - Single LYSO Gmsh meshed(Accepts Tile config -TileV0); 13 - 16 bars Gmsh meshed
    G4int GetStepSize() const {return StepSize;}		// Get Step size for particle propagation (user given)
    G4int GetSiPMmaterial() const {return SiPMmaterial;}	// Get value assignated to different possible materials for the SiPM
    G4int GetScintMat() const {return scint;}		// Get value assignated to different possible materials for the scintillator
	G4int GetESRTrue() const {return ESR;}			// Get value assignated to possible ESR cover (1 == True, 0 == False)
    G4int Getrndangle() const {return rndangle;}	// Get value assignated to possible random impact angles (1 == True, 0 == False)
    G4int GetYSymtrue() const {return NoYSym;}		// Get value assignated to possible Y axis symmetry in the Bar configuration (1 == True, 0 == False)
    G4int Getnrep() const {return nrep;}			// Get number of runs after the first one
    G4int GetZelem() const {return Zelem;}			// Get value assignated to number of element division across half a tesselated Bar
    G4int GetReflSiPM() const {return reflSiPM;}	// Get value assignated to possible Y axis symmetry in the Bar configuration (1 == True, 0 == False)
    G4int GetRnd_true() const {return RndGen[0];}	// If == 1 All random properties 
    G4int GetRnd_Part() const {return RndGen[1];}	// If == 1 Particle position set to random 
    G4int GetRnd_Geom() const {return RndGen[2];}	// If == 1 All geometry tolerances set to random
    G4int GetSpline() const {return Spline;}
    G4int GetFR4Refl() const {return FR4refl;}		// If == 1 the FR4 material has a fully reflective surface contact to the resin layer
    G4int GetforceBottomLine() const {return forceBottomLine;}	// REQUIRED FOR TILE CONFIG -TileV0
    G4int GetTile() const {return Tile;}			// If == 1 create Tile (Requires -GeomConfig 11 -TileV0 -ForceBottomLine)
    G4int GetnResinMach() const {return ResinMach;}	// If == 1 Machines a gap btw sipms
    G4int GetESRFinish() const {return ESRFinish;}	// If == 
    void matchSiPMf();								// If == 1 sets the width of the SiPM equal to the width of the LYSO
    G4int GetZnode() const {return Znode;}			// Defines the number of divisions in Z for GC 11/13 and the Gmsh meshed elements
    G4int GetnX() const {return nX;}				// Defines number of nodes along X for the -TileV0 configuration (default set to )
    G4int Getnodesec() const {return nodesec;}		// Defines the number of sections along Y for the -TileV0 configuration (default set to )
    G4int GetYstr() const {return Ystr;}			// Returns the value for each input node for the LYSO creation
    G4int GetGmshView() const {return GmshView;}	// If == 1 forces a gmsh visualization rather than a G4 visualization of the LYSO geom
    G4int GetOnode() const {return Onode;}

   
    // Double construction property Getters: geometry dimensions
    G4double GetGeom_DET_L() const {return DET_L;}	// Half length of the SiPM DETector thickness along Z in the bar configuration. Cmd -> 
    G4double GetGeom_RESIN_L() const {return RESIN_Z;}	// Half length of the SiPM resin package thickness along Z (without counting the SiPM thickness in Z) Cmd -> 
    G4double GetGeom_Resin_width() const {return Geom_Resin[1];}	// Half length of the SiPM resin package thickness along Y (direction of the SiPM package) centered in the middle bar Cmd -> 
    G4double GetGeom_LYSO_L() const {return Geom_LYSO[2];}	// Half length of the LYSO along Z Cmd -> 
    G4double GetGeom_LYSO_thick() const {return Geom_LYSO[1];}	// Half length of the LYSO along Y for the bar configuration Cmd -> 
    G4double GetGeom_DET_T() const {return DET_T;} // Half length of the detector SiPM Z thickness
    G4double GetGeom_DET_TX() const {return DET_TX;} //  Half length of the detector SiPM X width
    G4double GetGeom_DET_TX_tol() const {return DET_TX_tol;} //  X tolerance for the SiPM positioning, requires random geometry -rnd 0 1 0 or -rnd 1
    G4double GetGeom_DET_TY_tol() const {return DET_TY_tol;} //  Y tolerance for the SiPM positioning, requires random geometry -rnd 0 1 0 or -rnd 1
    G4double GetGeom_RESIN_H() const {return RESIN_H;} //  Half length of the resin Z width
    G4double GetGeom_RESIN_Y() const {return RESIN_Y;} //  Half length of the resin Y width
    G4double GetGeom_SiPM_Y() const {return SiPM_Y;} //  Half length of the detector SiPM height in Y
    G4double Get_GLUE_Y() const {return Glue_Y;} //  Half length of the detector glue y height 
    G4double GetXvec(int ind) const {return xv[ind];}
    G4double GetYvec(int ind) const {return yv[ind];}
    G4double GetYvecincr(int ind) const {return yvincr[ind];}
    G4double* GetYincr() const {return yincr;}
    G4double GetSZloc(){return SZ_loc;} // TileV0, location of the SiPm along the bottom part of the LYSO crystal from 0 to 1 as a percentage of Zhalf length
    G4double GetGlueZ() const {return Glue_Z;}

    
    // Int operational Getters : related to calculation of properties to return and/or variables to return
    G4int GetNSGAII() const {return NSGAII;} // If == 1 we prepare the NSGAII optimization
    G4int GetRootCreate() const {return RootCreate;}
    G4int GetTree_Hits() const {return MainTrees[0];}	//  If == 1  we store all the information about the photons hitting the SiPMs
    G4int GetTree_Detected() const {return MainTrees[1];} //  If == 1  we store all the information about the photons detected by the SiPMs
    G4int GetTree_Stepping() const {return MainTrees[2];} //  If == 1  we store all the information about all the photon steps 
    G4int GetTree_Tracking() const {return MainTrees[3];} //  If == 1  we store all the information about the photons tracks
    G4int GetTree_EndOfEvent() const {return MainTrees[4];} //  If == 1  we store all the information about the end of events
    G4int GetTree_EndOfRun() const {return MainTrees[5];} //  If == 1  we store all the information about the end of the run
    G4int SaveMesh() const {return SaveSTL;} //  If == 1  save the LYSO mesh as an STL file
    G4int GetVis() const {return VisTrue;} //  If == 1 , and visualization activated!!! we can see the G4 visualization
    G4int GetTimeTrue() const {return TimeTrue;}  //  If == 1  we 
    G4int GetKillTLTrue() const {return KillLTTrue;}
    G4int GetnEvents() const {return nEvents;}
    G4int GetLO() const {return ArgLO;} //  Returns the LO values
    G4int GetCT() const {return ArgCrossTalk;}  //  Returns the Cross Talk values
    G4int GetTP() const {return TotPh;} 
    G4int GetNPhotL() const {return nPhotL;} //  Returns photon hits in the Z<0 SiPM
    G4int GetNPhotR() const {return nPhotR;}//  Returns photon hits in the Z>0 SiPM
    G4int GetPhHits() const {return PhHit;}//  Returns total photon hits
    G4int GetNEdep() const {return NEdep;} //  Returns photon Energy deposition
    int GetRunEvt() const {return runevt;} //  Returns number of events per run
    G4int GetMuonFlag() const {return Muon;} //  if == 1 activates muon impacts rather than gamma 0.511 particles

    // Double operational Getters : related to calculation of properties to return
    G4double GetLYSO_Yield() const {return LYSOProps[0];} // Define material property: LYSO yield
    G4double GetLYSO_ScaleResolution() const {return LYSOProps[1];} // Define material property: LYSO scale resolution
    G4double GetLYSO_RiseT() const {return LYSOProps[2];} // Define material property: LYSO rise time 1
    G4double GetLYSO_DecayT() const {return LYSOProps[3];}// Define material property: LYSO decay time 1
    G4double GetKillTL() const {return KillLTime;} // Define event property: absolute time to kill all photons (by default deactivated)
    G4double GetNPhotTiming() const {return NPhotTiming;} 
    G4double GetPhotTiming();
    G4double GetGunX(G4int evt) const {return nGunPosX[evt];} // Define event property: X location of the impacting particle starting point (by default 0m)
    G4double GetGunY(G4int evt) const {return nGunPosY[evt];} // Define event property: Y location of the impacting particle starting point (byt default 0.5m)
    G4double GetnResinMachN() const {return ResinMachN;} // Define geometry property: do we machine a gap btw SiPMs? 1 = yes
    G4double GetMuonLYSOTrackLength() const { return MuonLYSOTrackLength;} // Get path length of the muon impacting particle
    G4double GetEdep() const {return Edep;} // Get Energy DEPosition of the impacting particle in the LYSO
    G4double GetMuonEdep() const {return MuonEdep;} // Get Energy DEPosition of the MUON impacting particle in the LYSO
    G4double GetPartDir(G4int i){return PartDir[i];} // Get Initial direction of the impacting particle along axis [i] (by default [0 -1 0])
    G4double GetPartXDispl(){return PartDisplX;} // Define event property: X displacement of the initial location of the impacting particle in case it has an angled impact
    G4double GetEvtLD(G4int evt) const {return nEventLD[evt];}   // Get Light Deposition (LD / Light Collection (LC)) in the SiPMs (total number of photons impacting and detected)
    G4double GetEvtLO(G4int evt) const {return nEventLO[evt];}  // Get Light Output(LO) in the SiPMs in the event  [evt]
    G4double GetEvtLSt(G4int evt) const {return nEventLSt[evt];}  // Get Light Output per track length of the impacting particle (LST) in the LYSO in the event  [evt]
    G4double GetEvtTim(G4int evt) const {return nEventLO[evt];}  // Get estimated timing for the event [evt]
    G4double GetLOIQR() const {return IQRLO;}  // Get LO IQR btw all events in the run
    G4double GetLDIQR() const {return IQRLD;}  // Get LD IQR btw all events in the run
    G4double GetLStIQR() const {return IQRLSt;}  // Get LST IQR btw all events in the run
    G4double GetLOAvg(G4int runid) const {return nRuntLOAvg[runid];}  // Get LO Average value btw all events in the run [runid]
    G4double GetLOStd(G4int runid) const {return nRuntLOStd[runid];}  // Get LO Standard Deviation value btw all events in the run [runid]
    G4double GetLOP50(G4int runid) const {return nRuntLOP50[runid];}  // Get LO P50 value btw all events in the run [runid]
    G4double GetLCP50(G4int runid) const {return nRuntLCP50[runid];}  // Get LD-LC Average value btw all events in the run [runid]
    G4double GetLStAvg(G4int runid) const {return nRuntLStAvg[runid];}  // Get LST Average value btw all events in the run [runid]
    G4double GetLStP50(G4int runid) const {return nRuntLStP50[runid];}  // Get LST P50 btw all events in the run [runid]
    G4double GetLDAvg(G4int runid) const {return nRuntLDAvg[runid];}  // Get LD-LC Average value btw all events in the run [runid]
    G4double GetLDStd(G4int runid) const {return nRuntLDStd[runid];}  // Get LD-LC Standard deviation btw all events in the run [runid]
    G4double GetTimAvg(G4int runid) const {return nRunTimingAvg[runid];}  // Get timing Average value btw all events in the run [runid]
    G4double GetTimStd(G4int runid) const {return nRunTimingStd[runid];}  // Get timing Standard deviation btw all events in the run [runid]
    G4double GetnEvtEdep(G4int runid) const {return nEdepEvts[runid];}  // 
    G4double GetVolume() const {return Volume;} // Get LYSO total volume
    
    // Double space exploration Getters: helped in the design exploration of the bar configurations
    G4double GetIncr() const {return incr;} // Test: [incr] forces an increase of thickness in the middle of the crystal equal to incr*3mm
    G4double GetIncrS() const {return incrS;} // Test: [incrS] forces an increase of thickness in the edge of the crystal equal to incr*3mm
    G4double GetIncrV() const {return incrV;} // Test: [incrV] forces an increase of thickness in the middle of the crystal equal to incr*3mm and equal decrease in the edge. Constant volume.
    G4double Getrad2Y() const {return rad2Y;}
    G4double Get_TileScale() const {return Tile_Scale;} // -TileV0 variation of the geometry based on a linear slope given by the scale from the nodes in the middle of the crystal 

    
	//*** Operational functions ***//
	// Adders: add value to a given stored value in the G4Args class
    void AddNEdep(){NEdep += 1;} // Add +1 to the Energy deposition events for the impacting particle
    void AddMuonEdep(G4double Edepadd){MuonEdep += Edepadd;}	// Add [Edepadd] to the event Muon Edep
    void AddEdep(G4double Edepadd){Edep += Edepadd;} // Add [Edepadd] to the event Edep
    void AddMuonLYSOTrackLength(G4double TLadd){MuonLYSOTrackLength += TLadd;} // Add [TLadd] to the impacting particle track length in the LYSO
    void AddPhHit(){PhHit += 1;} // +1 photon hit in a SiPM
    void AddLO(){ArgLO += 1;} // +1 photon hit in a SiPM for the LO calculation
    void AddCT(){ArgCrossTalk += 1;} // +1 photon hit in a SiPM not in the impacted LYSO edges (Cross-talk) - Requires GC 3/13
    void AddTP(){TotPh += 1;}
    void AddPhotR(){nPhotR += 1;}	// +1 photon hit in a SiPM in positive Z location
    void AddPhotL(){nPhotL += 1;}    // +1 photon hit in a SiPM in negative Z location
    void AddPhotTiming(G4double , G4double);
    void AddVolume(G4double Vol){Volume += Vol;}// LYSO volume addition operation for the Gmsh mesh volume calculation
    
	// Fillers: fill a list with each of the values in each 'event of a run' / 'run' 
    void FillAvgTim(G4int); // fill Average timing per run
    void FillAvgLO(G4int);  // fill Average LO per run
    void FillStdTim(G4int); // fill P50 timing per run
    void FillStdLO(G4int);  // fill P50 LO per run
    void FillEvtTim(G4int evt, G4double val){nEventTiming[evt]=val;}  // Arrival time storage
    void FillEvtLO(G4int evt, G4double val){nEventLO[evt]=val;}  // LO per event storage
    void FillEvtLSt(G4int evt, G4double val){nEventLSt[evt]=val;}  // LSP per event storage
    void FillEvtLD(G4int evt, G4double val){nEventLD[evt]=val;}  // LC per event storage
    
    //*** Setters: Used in follow up classes ***//
    void SetGeom_GLUE_Y(G4double val) { Glue_Y=val;}
    G4int SetVolume(G4double V){Volume=V;}

	//*** Currently under testing/Unused/to document ***//
    void GeomReinit();
    G4int FindEvents(G4String);
    void DefaultRadiusVect();
    void SetRadiusVect(G4double*, G4int, G4int);
    G4double GetGeomIndv(G4int runid) const {return RndGenIndv[runid];}  
	bool IsVolumeInList(const G4LogicalVolume* volume);
    void PushfScoringVolumeVec(G4LogicalVolume *LYSOTet_Logic) {fScoringVolumeVec.push_back(LYSOTet_Logic);}  
    G4double* GetNodeRadValues() const {return xv0;}

    
private:

    // Default values modifiable by arguments and able to be returned!!!
    G4int Oin=0, Rin=0, nrep=0; 
    G4int MainTrees[6]={0, 1, 0, 0, 1, 1};// Options to write(1)/orNot(0) the different output trees {Arrivals,Detected,Stepping,Tracking,EndOfEvent}
    G4double LYSOProps[4]={40000.,0.,60.,39.1};//Options to modify default LYSO properties {Yield,ScaleResolution,RiseTime,DecayTime}
    G4int RndGen[3]={1,1,1};             // Options regarding the random generator {Init,Particle,Geometry}
    G4double Geom_LYSO[3]={3./2.,3./2.,57./2.};//Options to modify default Geom properties {...}
    G4double Geom_Resin[3]={3./2.,51.5/2.,57./2.};//Options to modify default Geom properties {...}
    G4String OutName, MacName;
    G4double Detection[2]={3.5,1};//Options to modify photon detection
    G4String DefOutName="DefaultOutputName_Run"; // Default name for root files
    G4String DefRootFolder="./Results/"; // Default folder to store root files
    G4double KillSim[4]={0.,0.,0.,0.}; // Options on how to kill the simulation {Method(0 = all tracks have been killed, 1= nphotons have been detected at each SiPM, 2= x time has passed, 3= kills all photons after x local time, 4= use all parameters to kill the simulation) ,nphotons, global time, local time...}
    G4double StepSize=0.;    
    G4int GeomConfig=11; // Default geometric configuuration, [1, 2, 3, 11, 13] available, see top
    G4int VisTrue=1;
    G4int TimeTrue=1;
    G4int NPhotTiming=10; 
    G4double AvgTiming; // Storage of average photon timing
    G4double PhotTiming[2]={10.,0.};
    G4double TListR[10],TListL[10];
    G4double KillLTime=200; // time to kill all photons if activated
    G4int KillLTTrue=0; // kill all photons at killtime? 1 true
    G4int nEvents=0; // number of envets run along the simulation
    G4double *nEventTiming, *nEventLO,*nEventLSt, *nEventLD, *nRunTimingAvg, *nRuntLOAvg,*nRuntLOP50,*nRuntLCP50,*nRuntLStAvg,*nRuntLStP50, *nRuntLDAvg, *nRunTimingStd, *nRuntLOStd, *nRuntLDStd, *nGunPosX, *nGunPosY, *xv, *yv, *yincr, *yvincr;        
    G4int *nEdepEvts; // number of events with energy deposition
    G4double incr=0; // special study case, see top
    G4double incrS=0; // special study case, see top
    G4double incrV=0; // special study case, see top
    G4int Znode=1; // number of nodes for geoms 11 and 13
    G4int nX=3; // number of nodes in X for -TilveV0
    G4int nodesec=4; // number of sections for -TileV0
    G4int Onode=5;
    G4double* xv0 = NULL;
    G4int Muon = 0; // do we use Muons? 1 = True
    G4int dateflag;
    G4int RootCreate=1;
    G4int NSGAII=0; 
    G4int rad2Y=0;

    G4int ArgLO=0;
    G4int ArgCrossTalk=0;
    G4int TotPh=0; 
    G4int PhHit=0;
    G4double Edep=0.; // total energy deposition
    G4double MuonEdep=0.;
    G4int nPhotR=0; // number of photons in the +Z SiPM
    G4int nPhotL=0; // number of photons in the -Z SiPM
    G4int NEdep=0;
    G4int Ystr=0;
    G4double Volume=0.;
    G4int ResinMach=0; // machining holes between sipms to simulate single bar configuration (geometry modification)
    G4double ResinMachN=-1;
    G4int SiPMmatch=0; // match SipM width to the LYSO width?
    G4int SiPMmaterial=2;
    G4int AreaCte=0; // Change area with length of lyso automatically
    G4int Spline=1;
    G4int reflSiPM=0; 

    G4double DET_YMAX=3;
    G4double DET_XMAX=3;

    G4double DET_L=0.3/2.; // SiPM thickness default

    G4double DET_T=3./2; // SiPM height default
    G4double RESIN_Z=0.25; // resin thickness on thop of the SiPM by default

    G4double DET_TX=3./2; // SiPM widht in X by default
    G4double DET_TX_tol=0.;
    G4double DET_TY_tol=0.;
    G4double Ip=0.;

    G4double Glue_Y=3.; // Glue layer width by default
    G4double Glue_Z=0.2; // Glue layer thickness by default
    G4double RESIN_H=6.5/2; // Resin height by default
    G4double RESIN_Y=RESIN_H-0.5-Geom_LYSO[0]; // Actual location in Y of resin
    G4double SiPM_Y=+Geom_LYSO[0]+0.5-RESIN_H; // Actual SiPM Y location
    std::string YposStr;
    int runevt=0; // number of events to run
    char datechar [22];
    int Zelem=10; // mesh parameter along Z, number of elements
    int NoYSym=0; // for GeomConfig 11, no symmetry along the YX plane, both splines are independent and double of input values must be given in Ypos
    G4int FR4refl=0; // is the FR4 fully reflective?
    
    G4double PartDir[3];
    G4double PartDisplX=0; // Particle origin displacement in X at 0.5m in Y
    G4double PartAngle=0; // Angle of incidence of the particles
    G4int rndangle=0; // random impacting angle?
    G4int SaveSTL=0;
    G4int scint=1;
    G4int ESR=1; // Do we have ESR around the LYSO?
	G4int RndGenIndv[6] = {0};
	G4int forceBottomLine = 0; // needed for a flat bottom surface in -TileV0
	G4int Tile = 0;
    G4double SZ_loc=0.5; // location of the SiPM in TileV0 along the bottom surface
    G4int GmshView = 0; // load the Gmsh visualization
    G4double MuonLYSOTrackLength=0;
    G4int ESRFinish=0;
	std::vector<G4LogicalVolume*> fScoringVolumeVec;
    G4double IQRLO = 0; // IQR for LO
    G4double IQRLD = 0; // IQR for LC
    G4double IQRLSt = 0; // IQR for LSP
    G4int gunmesh=0; // do we use a uniform pattern of impacts?
    G4double Tile_Scale=0;
    G4String rootfolder;
    int imax=0, jmax=0; // number of impacts for the uniform impact pattern along X and Y

};    

#endif
