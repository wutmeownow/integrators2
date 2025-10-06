#include "ndcrescent_plotting.h"

int main(int argc, char* argv[]) {
    // Plot results for 10D, 5D, and 3D hyperspheres
    // Check if the user provided the correct number of arguments
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << "<double> <double> <double>\n";
        return 1; // Exit with error code
    }
    // Parse the command-line arguments
    double r1 = std::atof(argv[1]); // radius of first sphere centered at 'a'/2 from the origin
    double r2 = std::atof(argv[2]); // radius of second sphere displaced 'a'/2 away from origin
    double a = std::atof(argv[3]); // offset of spheres

    std::cout << "(r1,r2,a): " << r1 << " " << r2 <<" " << a << std::endl;

    int N; // number of points to generate in each trial
    std::vector<double> n_points_sqrt(19); // will have 19 different numbers of points from 2^6 to 2^24

    // 10D hyperspheres
    std::vector<double> results_10D(3); // vector to hold results (p, volume, stdev)
    std::vector<double> volume_10D(19); // vector to hold volume results for plotting 10
    std::vector<double> vol_err_10D(19); // vector to hold volume results for plotting

     // 5D hyperspheres
    std::vector<double> results_5D(3); // vector to hold results (p, volume, stdev)
    std::vector<double> volume_5D(19); // vector to hold volume results for plotting 10
    std::vector<double> vol_err_5D(19); // vector to hold volume results for plotting

    // 3D spheres
    std::vector<double> results_3D(3); // vector to hold results (p, volume, stdev)
    std::vector<double> volume_3D(19); // vector to hold volume results for plotting 10
    std::vector<double> vol_err_3D(19); // vector to hold volume results for plotting

    for (int i=6; i<=24; i++) {
        N = std::pow(2,i);
        n_points_sqrt[i-6] = std::sqrt(N);

        // 10D
        results_10D = monte_volume(10, N, r1, r2, a);
        volume_10D[i-6] = results_10D[1];
        vol_err_10D[i-6] = results_10D[2];
        // std::cout<<"10D volume: "<<results_10D[1]<<"\n"<<std::endl;

        // 5D
        results_5D = monte_volume(5, N, r1, r2, a);
        volume_5D[i-6] = results_5D[1];
        vol_err_5D[i-6] = results_5D[2];

        // 3D
        results_3D = monte_volume(3, N, r1, r2, a);
        volume_3D[i-6] = results_3D[1];
        vol_err_3D[i-6] = results_3D[2];
    }

    std::vector<double> zeros(19,0.0);

    // A multi panel plot
    auto canvas = new TCanvas("canvas","Canvas");
    canvas->Divide(3,2);  // divide in to 2x3 panels, plot volume on top of error
    
    // 10D plots
    canvas->cd(1);
    auto graph_10D = new TGraphErrors(19,n_points_sqrt.data(),volume_10D.data(),zeros.data(),vol_err_10D.data());
    graph_10D->SetTitle("10D Spheres;sqrt(N);Volume Shared");
    graph_10D->Draw();
    // canvas->cd(1)->SetLogy();
    canvas->cd(1)->SetLogx();

    canvas->cd(4);
    auto graph_10D_err = new TGraph(19,n_points_sqrt.data(),vol_err_10D.data());
    graph_10D_err->SetTitle(";sqrt(N);Error");
    graph_10D_err->Draw();
    canvas->cd(4)->SetLogy();
    canvas->cd(4)->SetLogx();

    // 5D plots
    canvas->cd(2);
    auto graph_5D = new TGraphErrors(19,n_points_sqrt.data(),volume_5D.data(),zeros.data(),vol_err_5D.data());
    graph_5D->SetTitle("5D Spheres;sqrt(N);Volume Shared");
    graph_5D->Draw();
    // canvas->cd(2)->SetLogy();
    canvas->cd(2)->SetLogx();

    canvas->cd(5);
    auto graph_5D_err = new TGraph(19,n_points_sqrt.data(),vol_err_5D.data());
    graph_5D_err->SetTitle(";sqrt(N);Error");
    graph_5D_err->Draw();
    canvas->cd(5)->SetLogy();
    canvas->cd(5)->SetLogx();

    // 3D plots
    canvas->cd(3);
    auto graph_3D = new TGraphErrors(19,n_points_sqrt.data(),volume_3D.data(),zeros.data(),vol_err_3D.data());
    graph_3D->SetTitle("3D Spheres;sqrt(N);Volume Shared");
    graph_3D->Draw();
    canvas->cd(3)->SetLogx();

    canvas->cd(6);
    auto graph_3D_err = new TGraph(19,n_points_sqrt.data(),vol_err_3D.data());
    graph_3D_err->SetTitle(";sqrt(N);Error");
    graph_3D_err->Draw();
    canvas->cd(6)->SetLogy();
    canvas->cd(6)->SetLogx();

    canvas->SaveAs("convergence.png");
    return 0; // Exit successfully
}