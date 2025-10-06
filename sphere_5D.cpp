#include "sphere_5D.h"

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

bool check_inside(int dim, std::vector<double>& v, double r) {
    // check if vector is inside a dim dimensional sphere centered at the origin
    std::vector<double> temp(v); // create a deep copy of the vector
    double mag = vector_mag(dim, temp);
    return mag <= r;
}

void add_vec(int dim, std::vector<double>& v1, std::vector<double>& v2) {
    // add components of vec2 to vec1 in place
    for (int j=0; j<dim; ++j){
        v1[j] = v1[j] + v2[j];
    }
}

int uni_loop(int dim, long int n_axis, int axis, double r, double incr, std::vector<double>& pos) {
    // axis is the current axis that is being incremented (ex: 0 is the first axis, dim-1 is the final axis)
    std::vector<double> temp(pos); // create a deep copy of the vector
    int sum = 0;
    if (axis<=(dim-1)) {
        // not outside the number of dimensions, call dim_loop again, modifying the current axis each time
        for (int i = 0; i<n_axis; i++){
            if (i==0) {
                sum += uni_loop(dim, n_axis, axis+1, r, incr, temp);
            }else {
                temp[axis] = temp[axis] + incr;
                sum += uni_loop(dim, n_axis, axis+1, r, incr, temp);
            }
        }
        return sum;
    } 

    // axis is not within number of dimensions, check if the given vector is inside d-sphere
    // helpers::print_double_vec(dim, pos);
    if (check_inside(dim, pos, r)) {
        return 1;
    }
    return 0;

}

std::vector<double> uniform_mc(int dim, int N, int r, long int n_axis) {
    // do mc integration of a dim dimensional sphere with radius r with N total points, n_axis points along each axis
    std::vector<double> results(2);
    std::vector<double> uniform_pos(dim,-r); // position vector for uniform distribution method, start at (-r,-r,-r,-r,-r)
    double box_length = 2.*r; // length of box to surround the 5d sphere
    std::cout<<"box length: "<<box_length<<std::endl;
    double axis_incr = box_length/(n_axis-1); // amount to increment along each axis to hit N points inside box
    std::cout<<"axis incr: "<<axis_incr<<std::endl;
    double box_volume = std::pow(box_length, 1.*dim);
    int N_inside= uni_loop(dim, n_axis, 0, r, axis_incr, uniform_pos);
    double volume = ((1.*N_inside)/N) * box_volume;
    results[0] = volume;
    results[1] = 0.0;
    return results;
}

int main(int argc, char* argv[]) {
    // use fixed grid, quasirandom sobel numbers, and pseudo random numbers to get volume of a 5D unit sphere
    // Check if the user provided the correct number of arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<int> <int>\n";
        return 1; // Exit with error code
    }
    // Parse the command-line arguments
    int lower_pow = std::atof(argv[1]); // lower power of 10 to start generating points at
    int upper_pow = std::atof(argv[2]); // upper power of 10 to generate points at

    std::cout << "(N_i, N_f): " << "10^"<<lower_pow<< " " <<"10^"<<upper_pow<< std::endl;

    int d = 5; // number of dimensions
    double r = 1.; // radius of 5D sphere
    long int N; // number of points to generate in each trial
    int n_trials = upper_pow-lower_pow+1; // total number of trials
    std::vector<double> n_points_sqrt(n_trials); // sqrt of all different numbers of points generated in each trial
    std::vector<double> n_points_sqrt_uni(n_trials); // sqrt of all different numbers of points generated in each trial for uniform points

    // uniform distribution method
    
    long int N_uniform; // total number of uniform points in each trial
    std::vector<double> volume_uni(n_trials); // volumes for uniform distribution
    std::vector<double> std_dev_uni(n_trials); // std dev for uniform distribution
    std::vector<double> results_uni(2); // vector to temporarily hold results for uniform mc


    // sobel quasirandom method

    for (int i=lower_pow; i<=upper_pow; i++) {
        N = std::pow(10,i);
        std::cout<<"N: "<<N<<std::endl;

        long int n_axis = std::pow(N*1.,1./5.); // rounded number of points along an axis for the uniform distribution
        N_uniform = std::pow(n_axis, 5);
        std::cout<<"N for uniform points: "<<N_uniform<<std::endl;

        n_points_sqrt[i-lower_pow] = std::sqrt(N); // sqrt of N for quasi and pseudo rng
        n_points_sqrt_uni[i-lower_pow] = std::sqrt(N_uniform); // sqrt of N for uniform distribution
        
        // uniform distr
        std::cout<<"num per axis: "<<n_axis<<std::endl;
        results_uni = uniform_mc(d, N_uniform, r, n_axis);

        std::cout<<"Uniform Volume: "<<results_uni[0]<<"\n"<<std::endl;
    }

    // // A multi panel plot
    // auto canvas = new TCanvas("canvas","Canvas");
    // canvas->Divide(3,2);  // divide in to 2x3 panels, plot volume on top of error
    
    // // 10D plots
    // canvas->cd(1);
    // auto graph_10D = new TGraph(19,n_points_sqrt.data(),volume_10D.data());
    // graph_10D->SetTitle("10D Spheres;sqrt(N);Volume Shared");
    // graph_10D->Draw();
    // // canvas->cd(1)->SetLogy();
    // canvas->cd(1)->SetLogx();

    // canvas->cd(4);
    // auto graph_10D_err = new TGraph(19,n_points_sqrt.data(),vol_err_10D.data());
    // graph_10D_err->SetTitle(";sqrt(N);Error");
    // graph_10D_err->Draw();
    // canvas->cd(4)->SetLogy();
    // canvas->cd(4)->SetLogx();

    // // 5D plots
    // canvas->cd(2);
    // auto graph_5D = new TGraph(19,n_points_sqrt.data(),volume_5D.data());
    // graph_5D->SetTitle("5D Spheres;sqrt(N);Volume Shared");
    // graph_5D->Draw();
    // // canvas->cd(2)->SetLogy();
    // canvas->cd(2)->SetLogx();

    // canvas->cd(5);
    // auto graph_5D_err = new TGraph(19,n_points_sqrt.data(),vol_err_5D.data());
    // graph_5D_err->SetTitle(";sqrt(N);Std Dev");
    // graph_5D_err->Draw();
    // canvas->cd(5)->SetLogy();
    // canvas->cd(5)->SetLogx();

    // // 3D plots
    // canvas->cd(3);
    // auto graph_3D = new TGraph(19,n_points_sqrt.data(),volume_3D.data());
    // graph_3D->SetTitle("3D Spheres;sqrt(N);Volume Shared");
    // graph_3D->Draw();
    // canvas->cd(3)->SetLogx();

    // canvas->cd(6);
    // auto graph_3D_err = new TGraph(19,n_points_sqrt.data(),vol_err_3D.data());
    // graph_3D_err->SetTitle(";sqrt(N);Std Dev");
    // graph_3D_err->Draw();
    // canvas->cd(6)->SetLogy();
    // canvas->cd(6)->SetLogx();

    // canvas->SaveAs("methods.png");
    return 0; // Exit successfully
}