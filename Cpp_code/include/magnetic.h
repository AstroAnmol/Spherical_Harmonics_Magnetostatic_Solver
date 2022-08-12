#ifndef MAGNETIC_H
#define MAGNETIC_H

#include <Eigen/Dense>

class magnetic{
    private:
        //constants
        double a, sep; //radius of particles (m) and separation between particles (m)
        double susc; // magnetic susceptibility
        double hmag; // magneitc field mag (A/m)
        double alpha; // angle between the magnetic field and particles (degress)
        double mu0=4*M_PI*1e-07;
        int L; // number of multipoles used
        //functions
        double nchoosek(int n, int k);
        double mag_field(double r, double theta, double phi);
        double d_lpmn(int n, int m, double x);
    public:
        magnetic(double arr[]); // array for magnetic force parameters [a, c, susc, Hmag, angle, L]
        
};
#endif