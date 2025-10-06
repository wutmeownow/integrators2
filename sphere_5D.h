#include <iostream>
#include <cstdlib> // For std::atoi and std::atof
#include "Math/QuasiRandom.h"
#include <cassert>
#include <cmath>
#include <vector>


void scale_vec(int dim, std::vector<double>& v, double scale);

double vector_mag(int dim, std::vector<double>& v);

bool check_inside(int dim, std::vector<double>& v, double r);

void add_vec(int dim, std::vector<double>& v1, std::vector<double>& v2);

int uni_loop(int dim, long int n_axis, int axis, double r, double incr, std::vector<double>& pos);

std::vector<double> uniform_mc(int dim, int N, int r, int n_axis);

// helper functions to print vectors
namespace helpers {

void print_double_vec(int dim, std::vector<double>& v) {
    // print out the contents of a dim dimensional vector
    for (int i = 0; i<dim; ++i) {
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}

void print_int_vec(int dim, std::vector<int>& v) {
    // print out the contents of a dim dimensional vector
    for (int i = 0; i<dim; ++i) {
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}

}