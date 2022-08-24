#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include "spherical_harmonics.h"
#include <cmath>
#include <chrono>
//define associate legendre functions for cos(theta)
double lpmn_cos(int m, int n, double theta){
    return std::pow(-1, m)*std::assoc_legendre(n, m, std::cos(theta));
}

double d_lpmn_cos(int m, int n, double theta){
    return ((m-n-1)*lpmn_cos(m,n+1,theta) + (n+1)*std::cos(theta)*lpmn_cos(m,n,theta))/(-std::sin(theta));
}

int main() {
    double a=1;
    double sep=2;
    double hmag=1;
    double alpha=0;
    double susc=1;
    double L=10;

    Eigen::VectorXd arr(6);
    arr<< a, sep, susc, hmag, alpha, L;


    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();

    spherical_harmonics trial(arr);
    std::cout<<trial.f<<std::endl;

    // Get ending timepoint
    auto stop = std::chrono::high_resolution_clock::now();
 
    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 
    std::cout << "Time taken by function: "
         << duration.count() << " microseconds" << std::endl;

    
}