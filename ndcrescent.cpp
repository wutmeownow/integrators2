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
    // scale a positive vector to have components between -scale/2 and scale/2
    for (int i = 0; i<dim; ++i) {
        *(vec+i) = *(vec+i) * scale - (scale/2);
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

bool check_inside(int dim, double vec[], double r, double offset) {
    // check if vector is inside a dim dimensional sphere offset from the origin in the first dimension
    vec[0] = vec[0] - offset;
    double offset_mag = vector_mag(dim, vec);
    vec[0] = vec[0] + offset;
    return offset_mag <= r;
}

double* monte_volume(int d, int N, double r1, double r2, double a) {
    QuasiRandomSobol srand(d); // sobol quasi random number generator in d dimensions
    double *d_rand = new double[d]; // double array to hold random vector in d dimensions
    double box_length = r1+r2+a; // length of box holding the two spheres
    std::cout<<"Box Length: "<<box_length<<std::endl;
    double box_volume = std::pow(box_length, d); // volume of the box

    int num_inside = 0;

    for (int i=0; i<N; ++i) {
        srand.RndmArray(d, d_rand); // generate a random vector
        // print_vec(d, d_rand);
        scale_vec(d, d_rand, box_length); // scale the vector by the box length
        // print_vec(d, d_rand);
        double d_rand_mag = vector_mag(d, d_rand); // calculate the magnitude of the vector
        if (check_inside(d, d_rand, r1, a/2) && check_inside(d, d_rand, r2, -a/2)) {
            // std::cout<<"Inside"<<std::endl;
            num_inside++;
        }
        // print_vec(d, d_rand);
        // std::cout<<d_rand_mag<<std::endl;
    }

    std::cout<<"n_inside: "<<num_inside<<std::endl;
    double p = num_inside/(N*1.);
    std::cout<<"p: "<<p<<std::endl;
    double volume = p*box_volume;
    std::cout<<"Volume: "<<volume<<std::endl;
    double stdev = 0;

    double* results = new double[3];
    results[0] = p;
    results[1] = volume;
    results[2] = stdev;

    delete[] d_rand;
    return results;
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

    double* results = monte_volume(d, N, r1, r2, a);
    double volume = results[1];
    double stdev = results[2];

    // *******

    // do not change the format below
    std::cout << "(r1,r2): " << r1 << " " << r2 << std::endl;
    std::cout << "(d,N,a): " << d << " " << N << " " << a << std::endl;
    std::cout << "volume: " << volume << std::endl;
    std::cout << "stat uncertainty: " << stdev << std::endl;

    delete[] results;
    return 0; // Exit successfully
}
