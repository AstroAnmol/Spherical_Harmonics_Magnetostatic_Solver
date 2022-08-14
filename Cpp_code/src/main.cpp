#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include "magnetic.h"
#include <cmath>
// define associate legendre functions for cos(theta)

int main() {
    double a=1;
    double sep=2;
    double hmag=1;
    double alpha=0;
    double susc=1;
    double L=10;

    Eigen::VectorXd arr(6);
    arr<< a, sep, susc, hmag, alpha, L;

    magnetic trial(arr);


    std::cout<< trial.f << std::endl<< std::endl;
    // std::cout<< trial.Beta1_1 << std::endl<< std::endl;
    // std::cout<< trial.Beta2_0 << std::endl<< std::endl;
    // std::cout<< trial.Beta2_1 << std::endl<< std::endl;

}