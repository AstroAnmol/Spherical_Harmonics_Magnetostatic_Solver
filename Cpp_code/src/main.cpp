#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include "magnetic.h"
#include <cmath>
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

    magnetic trial(arr);
    std::cout<<trial.f<<std::endl;

    // std::cout<< trial.mag_field(a, M_PI/3, M_PI/3).transpose()<< std::endl<< std::endl;
    // std::cout<< trial.mag_field(a, M_PI/5, M_PI/3).transpose()<< std::endl<< std::endl;
    // std::cout<< trial.mag_field(a, M_PI/3, M_PI/2).transpose()<< std::endl<< std::endl;
    // std::cout<< trial.mag_field(a, M_PI/8, M_PI/3).transpose()<< std::endl<< std::endl;
    
    // std::cout<<trial.integrand(M_PI/8, M_PI/7)<<std::endl;
    // std::cout<<trial.integrand(M_PI, M_PI/7)<<std::endl;
    // std::cout<<trial.integrand(3*M_PI/5, M_PI/9)<<std::endl;
    // std::cout<<trial.integrand(4*M_PI/9, M_PI/4)<<std::endl;

    
}