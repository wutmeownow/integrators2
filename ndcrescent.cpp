#include "ndcrescent.h"

using namespace ROOT::Math;

void print_vec(int dim, std::vector<double>& v) {
    // print out the contents of a dim dimensional vector
    for (int i = 0; i<dim; ++i) {
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}

void scale_vec(int dim, std::vector<double>& v, double scale) {
    // scale a positive vector to have components between -scale/2 and scale/2
    for (int i = 0; i<dim; ++i) {
        v[i] = v[i] * scale - (scale/2);
    }
}

double vector_mag(int dim, std::vector<double>& v) {
    // calculate the magnitude of a vector
    double sum = 0;
    for (int j=0; j<dim; ++j){
        sum += v[j]*v[j];
    }
    double mag = std::pow(sum, 0.5);
    return mag;   

}

bool check_inside(int dim, std::vector<double>& v, double r, double offset) {
    // check if vector is inside a dim dimensional sphere offset from the origin in the first dimension
    v[0] = v[0] - offset;
    double offset_mag = vector_mag(dim, v);
    v[0] = v[0] + offset;
    return offset_mag <= r;
}

void add_vec(int dim, std::vector<double>& v1, std::vector<double>& v2) {
    // add components of vec2 to vec1 in place
    for (int j=0; j<dim; ++j){
        v1[j] = v1[j] + v2[j];
    }
}


std::vector<double> monte_volume(int d, int N, double r1, double r2, double a) {
    std::cout << "Dimensionality (d): " << d << std::endl;
    assert(d > 0 && d < 10000);  // Add a sanity check


    QuasiRandomSobol qrsg(d); // sobol quasi random number generator in d dimensions
    std::vector<double> rand_pos(d); // double vector to hold random vector in d dimensions
    double box_length = r1+r2+a; // length of box holding the two spheres
    std::cout<<"Box Length: "<<box_length<<std::endl;
    double box_volume = std::pow(box_length, 1.*d); // volume of the box
    std::cout << "Box Volume: " << box_volume << std::endl;
    int num_inside = 0; // count how many random points land inside both spheres

    double sum_magsq = 0.0; // sum of the magnitude squared of all randomly generated vectors
    std::vector<double> sum_pos(d, 0.0); // vector which is the sum of all random vector components
    print_vec(d, sum_pos);

    double rand_mag;
    // go through N samples of the volume
    for (int i=0; i<N; ++i) {
        qrsg.RndmArray(1, rand_pos.data()); // generate a random vector
        // print_vec(d, rand_pos);
        scale_vec(d, rand_pos, box_length); // scale the vector by the box length
        // print_vec(d, rand_pos);

        rand_mag = vector_mag(d, rand_pos); // calculate the magnitude of the vector

        // add magnitude squared of this vector to d_rand_mag_sum
        sum_magsq += std::pow(rand_mag, 2.0);

        // add d_rand components to d_rand_sum array
        add_vec(d, sum_pos, rand_pos);
        // print_vec(d, d_rand_sum);
        
        if (check_inside(d, rand_pos, r1, a/2) && check_inside(d, rand_pos, r2, -a/2)) {
            // std::cout<<"Inside"<<std::endl;
            num_inside++;
        }
        // print_vec(d, d_rand);
        // std::cout<<d_rand_mag<<std::endl;
    }

    // calculate average vector
    for (int j=0; j<d; ++j){
        sum_pos[j] = sum_pos[j]/(1.* N);
    }
    double exp_xsqr = std::pow(vector_mag(d, sum_pos), 2);

    std::cout<<"n_inside: "<<num_inside<<std::endl;
    double p = num_inside/(N*1.);
    std::cout<<"p: "<<p<<std::endl;
    double volume = p*box_volume;
    std::cout<<"Volume: "<<volume<<std::endl;
    std::cout<<"mag sqr sum "<<sum_magsq<<std::endl;

    double stdev = std::pow((sum_magsq/(1.*N) - exp_xsqr), 0.5) / std::pow(N, 0.5);
    std::cout<<"std dev: "<<stdev<<std::endl;
    std::vector<double> results(3);
    results[0] = p;
    results[1] = volume;
    results[2] = stdev;

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
    std::cout << "Dimensionality (d): " << d << std::endl;
    assert(d > 0 && d < 10000);  // Add a sanity check

    std::vector<double> results = monte_volume(d, N, r1, r2, a);
    double volume = results[1];
    double stdev = results[2];

    // *******

    // do not change the format below
    std::cout << "(r1,r2): " << r1 << " " << r2 << std::endl;
    std::cout << "(d,N,a): " << d << " " << N << " " << a << std::endl;
    std::cout << "volume: " << volume << std::endl;
    std::cout << "stat uncertainty: " << stdev << std::endl;

    return 0; // Exit successfully
}
