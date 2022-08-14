#include <iostream>
#include <cmath>
#include "magnetic.h"
#include <eigen-3.4.0/Eigen/Dense>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__
// Intiator
magnetic::magnetic(Eigen::VectorXd arr){
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

    H0<< H_perp, 0, H_prll;

    for (int m= 0; m < 2; m++){
        Eigen::MatrixXd X(L,L), Delta_m(L,L), Gamma_m(L,L); 
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
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
        Am.block(0,0,L,L)=X;
        Am.block(0,L,L,L)=Delta_m;
        Am.block(L,0,L,L)=Gamma_m;
        Am.block(L,L,L,L)=X;

        //qm vector
        Eigen::VectorXd qm(L);
        if (m==0){
            qm(0)=-H_prll*std::pow(a,3)*(1-mu/mu0);
        }
        else if (m==1){
            qm(0)=H_perp*std::pow(a,3)*(1-mu/mu0);
        }
        
        //2L Q vector
        Eigen::VectorXd Qm(2*L);
        Qm.block(0,0,L,1)=qm;
        Qm.block(L,0,L,1)=qm;

        //solve linear system
        Eigen::VectorXd Beta_m(2*L);

        Beta_m=Am.colPivHouseholderQr().solve(Qm);
        if (m==0){
            Beta1_0=Beta_m.block(0,0,L,1);
            Beta2_0=Beta_m.block(L,0,L,1);
        }
        else if (m==1){
            Beta1_1=Beta_m.block(0,0,L,1);
            Beta2_1=Beta_m.block(L,0,L,1);
        }
    };
    std::cout<< "Linear System Solved"<<std::endl;
    
    //Create a 3D spherical mesh
    int N =180;
    double dang= M_PI/N;
    Eigen::VectorXd inc= Eigen::VectorXd::LinSpaced(N,dang/2, M_PI + dang).transpose();
    Eigen::VectorXd az= Eigen::VectorXd::LinSpaced(2*N,dang/2, 2*M_PI + dang).transpose();

    f=Eigen::Vector3d::Zero();

    // std::cout<<"for loop started"<<std::endl;
    // Formulating the Maxwell Stress Tensor in Spherical Coordinates
    for (int i = 0; i < 2*N; i++){
        double p;
        if (i==0 or i==(2*N-1)){
            p=1;}
        else if (i%2==0){
            p=2;}
        else{ p=4;}
        for (int j = 0; j < N; j++){
            double q;
            if (j==0 or j==(N-1)){
                q=1;}
            else if (j%2==0){
                q=2;}
            else{ q=4;}
            double ph= az[i];
            double th= inc[j];
            f=f+ a*p*q*integrand(th, ph);
        }
    }
    f=f*dang*dang/9.0;

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

Eigen::Vector3d magnetic::mag_field(double r, double theta, double phi){
    double Hr=0, Hth=0, Hphi=0;
    for (int l = 1; l < L+1; l++){
        for (int m = 0; m < 2; m++){
            double Hrs=0, Hths=0, Hphis=0;
            for (int s = m; s < L+1; s++){
                double Psm=lpmn_cos(m, s, theta);
                double dPsm=d_lpmn_cos(m, s, theta);
                double ls_choose_sm=nchoosek(l+s, s+m);
                double r_pow_s1=std::pow(r, (s-1));
                double sep_pow=std::pow(sep, l+s+1);
                double minus_one_pow=std::pow(-1, s+m);
                double r_pow_times_Psm=r_pow_s1*Psm;

                // r component
                double additional=(minus_one_pow)*s*r_pow_times_Psm/(sep_pow);
                double Hrs=Hrs + additional*(ls_choose_sm);
                // theta component
                double Hths=Hths + (minus_one_pow)*ls_choose_sm*(r_pow_s1*dPsm)/(sep_pow);
                // phi component
                if (m==1){
                    Hphis= Hphis + (minus_one_pow)*ls_choose_sm*r_pow_times_Psm/(std::sin(theta))/(sep_pow);
                }

                double Plm=lpmn_cos(m, l, theta);
                double dPlm=d_lpmn_cos(m, l, theta);
                double r_pow_l2=std::pow(r, l+2);

                if (m==0){
                    // R component
                    Hr=Hr + (((l+1)*Beta1_0[l-1]*(Plm/r_pow_l2) -  Beta2_0[l-1]*Hrs)*std::cos(m*phi));
                    // Theta component
                    Hth=Hth + ((Beta1_0[l-1]*(dPlm/r_pow_l2) + Beta2_0[l-1]*Hths)*std::cos(m*phi));
                }
                else if (m==1){
                    // R Component
                    Hr=Hr + (((l+1)*Beta1_1[l-1]*(Plm/r_pow_l2) - Beta2_1[l-1]*Hrs)*std::cos(m*phi));
                    // Theta component
                    Hth=Hth + ((Beta1_1[l-1]*(dPlm/r_pow_l2) + Beta2_1[l-1]*Hths)*std::cos(m*phi));
                    // Phi component
                    Hphi=Hphi + ((Beta1_1[l-1]*(Plm/(std::sin(theta)/r_pow_l2)) + Beta2_1[l-1]*Hphis)*std::sin(phi));
                }
            }
        }
    }
    Hth=-Hth;
    Eigen::Vector3d magfield;
    magfield << Hr, Hth, Hphi;
    return magfield;
}

//integrand function
Eigen::Vector3d magnetic::integrand(double th, double ph){
    //transformation matrix
    Eigen::Matrix3d pre, post, T_cart;
    pre<<   std::sin(th)*std::cos(ph), std::cos(th)*std::cos(ph), -std::sin(ph),
            std::sin(th)*std::sin(ph), std::cos(th)*std::sin(ph),  std::cos(ph),
            std::cos(th), -std::sin(th),  0;
    post= pre.transpose();
    Eigen::Vector3d H0_sph, H_sph, H_cart, rn_hat;
    H0_sph=post*H0;
    H_sph= mag_field(a, th, ph) + H0_sph;
    H_cart=pre*H_sph;
    double h=H_cart.norm();
    T_cart=mu0*(H_cart*H_cart.transpose() - 0.5*(h*h)*Eigen::Matrix3d::Identity());
    rn_hat<<std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),std::cos(th);
    return std::sin(th)*T_cart*rn_hat;
}

//define associate legendre functions for cos(theta)
double magnetic::lpmn_cos(int m, int n, double theta){
    return std::assoc_legendre(n, m, std::cos(theta));
}

double magnetic::d_lpmn_cos(int m, int n, double theta){
    return ((m-n-1)*lpmn_cos(m,n+1,theta) + (n+1)*std::cos(theta)*lpmn_cos(m,n,theta))/(-std::sin(theta));
}