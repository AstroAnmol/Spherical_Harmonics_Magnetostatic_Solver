#ifndef MAGNETIC_H
#define MAGNETIC_H
#include <cmath>
#include <eigen-3.4.0/Eigen/Dense>

class magnetic{
    private:
    public:
        magnetic(Eigen::VectorXd arr); // array for magnetic force parameters [a, c, susc, Hmag, angle, L]
        
};
#endif