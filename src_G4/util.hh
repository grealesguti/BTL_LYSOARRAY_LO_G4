#ifndef UTIL_HH
#define UTIL_HH

// Macro to suppress unused variable warnings
#define UNUSED(x) (void)(x)

// Include standard C and C++ headers for file I/O, string manipulation, and mathematical operations
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <algorithm> // std::swap
#include <vector> // std::vector
#include <fstream>
#include <limits>
#include <unistd.h> // for access()
#include <stdlib.h> // for getenv()
#include <stdio.h> // for sprintf(), fopen()
#include <ctime>
#include <G4Types.hh> // Geant4 types
#include <set>
#include <gmsh.h> // Gmsh API

#include <sstream>
#include <iterator>

// Function prototypes for utility functions
FILE *open_file(const char *filename, const char *mode); // Open a file
char *find_file(const char *filename); // Find a file in the system
int read_tsv_file(const char *filename, double *energy, double *values, double xunit, double yunit); // Read a TSV file
int partition(double arr[], int start, int end); // Partition function for quicksort
void quickSort(double arr[], int start, int end); // Quicksort algorithm
double LYSOMeshVolume(G4double* xv0, G4int Onode, G4int Znode); // Calculate volume of LYSO mesh
double *Str2DChar(std::string strinput, G4int nn); // Convert string to 2D char array
G4double TetraVolume(G4double *x,G4double *y,G4double *z); // Calculate volume of a tetrahedron

// Functions for manipulating vectors of vectors
std::vector<std::vector<double>> appendReversedInitialVector(const std::vector<std::vector<double>>& inputVector); // Append reversed initial vector
std::vector<std::vector<double>> reverseSubvectors(const std::vector<std::vector<double>>& inputVector); // Reverse subvectors
std::vector<std::vector<double>> combineAllIndices(const std::vector<double>& vec1, const std::vector<double>& vec2); // Combine all indices
std::vector<std::vector<int>> combineAllIndicesint(const std::vector<int>& vec1, const std::vector<int>& vec2); // Combine all indices (integer version)
std::vector<std::vector<double>> calculatePoints(const std::vector<double>& Y, double Xin, double Zmin, double Zmax, double Yzero); // Calculate points
std::vector<int> createGmshPoints(const std::vector<std::vector<double>>& points, int taginit); // Create Gmsh points
std::vector<int> createGmshLines(const std::vector<int>& , int , const std::set<int>& ); // Create Gmsh lines
int createGmshSurface(const std::vector<int>& lineTags, int taginit); // Create Gmsh surface
std::vector<double> extractRange(const std::vector<double>& Y_all, int sec, int nodesec); // Extract range from vector
void printVectorOfVectors(const std::vector<std::vector<int>>& vec); // Print vector of vectors
void printVector(const std::vector<int>& vec); // Print vector of int
void printVectordouble(const std::vector<double>& vec); // Print vector of doubles
std::vector<int> getElementsAtIndex(const std::vector<std::vector<int>>& vec, int index, int a, int b); // Get elements at index
std::vector<double> generateEquispacedSegments(double Xmin, double Xmax, int numElements); // Generate equispaced segments
std::vector<double> convertArrayToVector(const double* arr, size_t size); // Convert array to vector
std::vector<int> getOddElementsAtIndex(const std::vector<std::vector<int>>& vec, int index, int a, int b); // Get odd elements at index
std::vector<int> getEvenElementsAtIndex(const std::vector<std::vector<int>>& vec, int index, int a, int b); // Get even elements at index
std::vector<int> getValuesAtIndices(const std::vector<int>& vec, const std::vector<int>& indices); // Get values at indices

#endif
