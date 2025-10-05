#include <iostream>
#include <cstdlib> // For std::atoi and std::atof
#include "Math/QuasiRandom.h"
#include <cassert>
#include <cmath>
#include <vector>

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