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
    std::vector<double> temp(v); // create a deep copy of the vector
    temp[0] = temp[0] - offset;
    double offset_mag = vector_mag(dim, temp);
    return offset_mag <= r;
}

void add_vec(int dim, std::vector<double>& v1, std::vector<double>& v2) {
    // add components of vec2 to vec1 in place
    for (int j=0; j<dim; ++j){
        v1[j] = v1[j] + v2[j];
    }
}


std::vector<double> monte_volume(int d, int N, double r1, double r2, double a) {
    // std::cout << "Dimensionality (d): " << d << std::endl;
    assert(d > 0 && d < 10000);  // Add a sanity check


    QuasiRandomSobol sobol_rng(d); // sobol quasi random number generator in d dimensions
    std::vector<double> rand_pos(d); // double vector to hold random vector in d dimensions
    double box_length = r1+r2+a; // length of box holding the two spheres
    // std::cout<<"Box Length: "<<box_length<<std::endl;
    double box_volume = std::pow(box_length, 1.*d); // volume of the box
    // std::cout << "Box Volume: " << box_volume << std::endl;
    int num_inside = 0; // count how many random points land inside both spheres
    double stdev; // standard deviation of the points inside spheres

    double sum_magsq = 0.0; // sum of the magnitude squared of all randomly generated vectors inside both spheres
    std::vector<double> sum_pos(d, 0.0); // vector which is the sum of all random vector components inside both spheres
    // print_vec(d, sum_pos);

    double rand_mag;
    // go through N samples of the volume
    for (int i=0; i<N; ++i) {
        sobol_rng.RndmArray(1, rand_pos.data()); // generate a random vector
        // print_vec(d, rand_pos);
        scale_vec(d, rand_pos, box_length); // scale the vector by the box length
        // print_vec(d, rand_pos);


        rand_mag = vector_mag(d, rand_pos); // calculate the magnitude of the vector


        if (check_inside(d, rand_pos, r1, a/2) && check_inside(d, rand_pos, r2, -a/2)) {
            // std::cout<<"Inside"<<std::endl;
            num_inside++;

                    // add magnitude squared of this vector to d_rand_mag_sum
            sum_magsq += std::pow(rand_mag, 2.0);

            // add d_rand components to d_rand_sum array
            add_vec(d, sum_pos, rand_pos);
            // print_vec(d, d_rand_sum);
        }
        // print_vec(d, d_rand);
        // std::cout<<d_rand_mag<<std::endl;
    }

    if (num_inside <= 1) {
        stdev = 0.0;
    } else {
        // calculate average vector
        for (int j=0; j<d; ++j){
            sum_pos[j] = sum_pos[j]/(1.* num_inside);
        }
        double exp_pos_sqr = std::pow(vector_mag(d, sum_pos), 2); // magnituded squared of average vector inside both spheres
        stdev = std::pow((sum_magsq/(1.*num_inside) - exp_pos_sqr), 0.5) / std::pow(num_inside-1, 0.5); // standard deviation in points generated inside both spheres
    }
    

    // std::cout<<"n_inside: "<<num_inside<<std::endl;
    double p = num_inside/(N*1.);
    // std::cout<<"p: "<<p<<std::endl;
    double volume = p*box_volume;
    // std::cout<<"Volume: "<<volume<<std::endl;
    // std::cout<<"\n"<<std::endl;
    // std::cout<<"mag sqr sum "<<sum_magsq<<std::endl;
    
    // stdev = stdev * (box_volume / N); // scale this standard error by the total box volume divided by total number of points generated
    // std::cout<<"std dev: "<<stdev<<std::endl;
    std::vector<double> results(3);
    results[0] = p;
    results[1] = volume;
    results[2] = stdev;

    return results;
}


