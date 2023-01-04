#ifndef MAG_FORCE_H
#define MAG_FORCE_H
#include <cmath>
#include <eigen-3.4.0/Eigen/Dense>

class mag_force{
    private:
        int N; // number of particles
        double susc; //susceptibility
        double a; //radius of particles
        double susc_eff=3*susc/(susc+3); // effective susceptiblity
        Eigen::Vector3d mag_field;
    // particle positions
        Eigen::VectorXd pos_x;
        Eigen::VectorXd pos_y;
        Eigen::VectorXd pos_z;
    // particle moments
        Eigen::VectorXd moment_x;
        Eigen::VectorXd moment_y;
        Eigen::VectorXd moment_z;

    // functions
        Eigen::Vector3d get_par_pos(int i);
        Eigen::Vector3d get_par_moment(int i);
        Eigen::Vector3d set_par_pos(int i, Eigen::Vector3d pos);
        Eigen::Vector3d set_par_moment(int i, Eigen::Vector3d mom);
    public:
        mag_force();
        void MDM();
        
};
#endif