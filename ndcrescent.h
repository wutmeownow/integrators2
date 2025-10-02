#include <iostream>
#include <cstdlib> // For std::atoi and std::atof
#include "TApplication.h"
#include "Math/QuasiRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGClient.h"
#include "TStyle.h"
#include <cmath>
#include <vector>

// print a dim dimensional vector
void print_vec(int dim, std::vector<double>& v);

// scale double vector in place
void scale_vec(int dim, std::vector<double>& v, double scale);

// calculate the magnitude of a double vector
double vector_mag(int dim, std::vector<double>& v);

// check if a vector defined off the origin is inside a sphere offset from the origin
bool check_inside(int dim, std::vector<double>& v, double r, double offset);

// add array vec2 to vec1 in place (modify vec1)
void add_vec(int dim, std::vector<double>& v1, std::vector<double>& v2);

// do monte carlo integration for two d-spheres - return p, the final volume, and std dev
std::vector<double> monte_volume(int d, int N, double r1, double r2, double a);