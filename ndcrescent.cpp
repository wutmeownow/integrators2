#include "ndcrescent.h"

using namespace ROOT::Math;

void print_vec(int dim, double vec[]) {
    // print out the contents of a dim dimensional vector
    for (int i = 0; i<dim; ++i) {
        std::cout<<vec[i]<<" ";
    }
    std::cout<<std::endl;
}

void scale_vec(int dim, double* vec, double scale) {
    // scale a vector
    for (int i = 0; i<dim; ++i) {
        *(vec+i) = *(vec+i) * scale;
    }
}

double vector_mag(int dim, double vec[]) {
    // calculate the magnitude of a vector
    double sum = 0;
    for (int j=0; j<dim; ++j){
        sum += vec[j]*vec[j];
    }
    double mag = std::pow(sum, 0.5);
    return mag;   

}

int main(int argc, char* argv[]) {
    // Check if the user provided the correct number of arguments
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <int1> <int2> <double> <double> <double>\n";
        return 1; // Exit with error code
    }

    // Parse the command-line arguments
    int d = std::atoi(argv[1]); // number of dimensions
    int N = std::atoi(argv[2]); // number of trials
    double r1 = std::atof(argv[3]); // radius of first sphere centered at 'a'/2 from the origin
    double r2 = std::atof(argv[4]); // radius of second sphere displaced 'a'/2 away from origin
    double a = std::atof(argv[5]); // offset of spheres
    
    // ******* add your code here 
    QuasiRandomSobol srand(d); // sobol quasi random number generator in d dimensions
    double d_rand[d]; // double array to hold random vector in d dimensions
    double box_length = r1+r2+a; // length of box holding the two spheres
    std::cout<<"Box Length: "<<box_length<<std::endl;
    double box_volume = std::pow(box_length, d); // volume of the box


    for (int i=0; i<N; ++i) {
        srand.RndmArray(2, d_rand); // generate a random vector
        print_vec(d, d_rand);
        scale_vec(d, d_rand, box_length); // scale the vector by the box length
        print_vec(d, d_rand);
        double d_rand_mag = vector_mag(d, d_rand); // calculate the magnitude of the vector
        std::cout<<d_rand_mag<<std::endl;
    }

    double volume = 0;
    double stdev = 0;

    // *******

    // do not change the format below
    std::cout << "(r1,r2): " << r1 << " " << r2 << std::endl;
    std::cout << "(d,N,a): " << d << " " << N << " " << a << std::endl;
    std::cout << "volume: " << volume << std::endl;
    std::cout << "stat uncertainty: " << stdev << std::endl;

    return 0; // Exit successfully
}
