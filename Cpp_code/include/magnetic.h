#ifndef MAGNETIC_H
#define MAGNETIC_H
#include <cmath>
#include <eigen-3.4.0/Eigen/Dense>

class magnetic{
    private:
        //constants
        double a, sep; //radius of particles (m) and separation between particles (m)
        double susc; // magnetic susceptibility
        double hmag; // magneitc field mag (A/m)
        double alpha; // angle between the magnetic field and particles (degress)
        double mu0=4*M_PI*1e-07;
        int L; // number of multipoles used
        Eigen::Vector3d H0;
        //functions
        double nchoosek(int n, int k);
        double lpmn_cos(int n, int m, double theta);
        double d_lpmn_cos(int n, int m, double theta);
    public:
        Eigen::Vector3d f;
        magnetic(Eigen::VectorXd arr); // array for magnetic force parameters [a, c, susc, Hmag, angle, L]
        Eigen::Vector3d integrand(double th, double ph);
        Eigen::Vector3d mag_field(double r, double theta, double phi);
        Eigen::VectorXd Beta1_0, Beta2_0, Beta1_1, Beta2_1;
};
#endif