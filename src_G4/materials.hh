// Include necessary Geant4 headers for material management
#include "G4NistManager.hh" // Class for managing NIST materials
#include "G4Material.hh" // Class representing a material in Geant4

// Function to create a custom material named "lyso" with specified light yield, rise time, and scale resolution
G4Material *get_lyso(double light_yield, double rise_time, double scale_resolution);

// Function to create a custom material named "BC400" with specified light yield, rise time, and scale resolution
G4Material *get_BC400(double light_yield, double rise_time, double scale_resolution);

// Function to create a custom material named "BC408" with specified light yield, rise time, and scale resolution
G4Material *get_BC408(double light_yield, double rise_time, double scale_resolution);

// Function to create a custom material named "rtv" without specifying light yield, rise time, and scale resolution
// This function likely uses default values or predefined properties for the "rtv" material
G4Material *get_rtv(void);
