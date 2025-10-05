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

    int N; // number of points to generate in each trial
    int n_trials = upper_pow-lower_pow+1; // total number of trials
    std::vector<double> n_points_sqrt(n_trials); // sqrt of all different numbers of points generated in each trial

    // sobel quasirandom method

    for (int i=lower_pow; i<=upper_pow; i++) {
        N = std::pow(10,i);
        n_points_sqrt[i-lower_pow] = std::sqrt(N);
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