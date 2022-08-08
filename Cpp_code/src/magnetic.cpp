#include <iostream>
#include <cmath>
#include "magnetic.h"
#include <Eigen/Dense>

// Intiator
magnetic::magnetic(double arr[]){
    // variable assignment
    a=arr[0];
    sep=arr[1]*a;
    susc=arr[2];
    hmag=arr[3];
    alpha=arr[4]*M_PI/180.0;
    L=arr[5];

    double mu=(1+susc)*mu0;

    double H_perp=hmag*std::sin(alpha);
    double H_prll=hmag*std::cos(alpha);

    

    for (int m= 0; m < 2; m++)
    {
        Eigen::MatrixXd X(L,L), Delta_m(L,L), Gamma_m(L,L); 
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                // X matrix
                if (i==j)
                {
                    X(i,j)=(i+1)*(mu/mu0) + (i+1) + 1;
                }
                // Delta and Gamma matrix
                Delta_m(i,j)=std::pow((-1),((i+1)+m))*((i+1)*(mu/mu0)-(i+1))*nchoosek(i+1+j+1, j+1-m)*std::pow(a,(2*(i+1)+1))/std::pow(sep,(i+1+j+1+1));
                Gamma_m(i,j)=std::pow((-1), (i+1+j+1))*Delta_m(i,j);
            }
            
        }
        // 2L X 2L Matrix
        Eigen::MatrixXd Am(2*L, 2*L);

        
    }
    
    
}


double magnetic::nchoosek(int n, int k){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}