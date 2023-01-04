#include <iostream>
#include <cmath>
#include "mag_force.h"
#include <eigen-3.4.0/Eigen/Dense>

mag_force::mag_force(){
    double sep=2.2;
    a=1;
    N=2;

    Eigen::VectorXd placeholder(2);
    placeholder = Eigen::VectorXd::Zero();

    pos_x=placeholder;
    pos_y=placeholder;
    pos_z=placeholder;
    
    moment_x=placeholder;
    moment_y=placeholder;
    moment_z=placeholder;

    Eigen::Vector3d pos_0(0,0,0);
    Eigen::Vector3d pos_1(0,0,sep);

    set_par_pos(0,pos_0);
    set_par_pos(1,pos_1);
    
}

void mag_force::MDM(){
    Eigen::Vector3d fixed_dipole=(4/3)*M_PI*(a*a*a)*susc_eff*mag_field;
    
    set_par_moment(0,fixed_dipole);
    set_par_moment(1,fixed_dipole);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Eigen::Vector3d C_ij= get_par_pos(i)-get_par_pos(j);
            double c_ij=C_ij.norm();

            Eigen::Vector3d H_dip;
            H_dip=1/(4*M_PI)*(3*C_ij*C_ij.dot(get_par_moment(i))/(std::pow(c_ij,5)) 
            - get_par_moment(i)/(std::pow(c_ij,3)));

            
        }
    }
    
}

Eigen::Vector3d mag_force::get_par_pos(int i){
    Eigen::Vector3d pos;
    pos[0]=pos_x[i];
    pos[1]=pos_y[i];
    pos[2]=pos_z[i];
    return pos;
}

Eigen::Vector3d mag_force::get_par_moment(int i){
    Eigen::Vector3d mom;
    mom[0]=moment_x[i];
    mom[1]=moment_y[i];
    mom[2]=moment_z[i];
    return mom;
}

Eigen::Vector3d mag_force::set_par_pos(int i, Eigen::Vector3d pos){
    pos_x[i]=pos[0];
    pos_y[i]=pos[1];
    pos_z[i]=pos[2];
}


Eigen::Vector3d mag_force::set_par_moment(int i, Eigen::Vector3d mom){
    moment_x[i]+=mom[0];
    moment_y[i]+=mom[1];
    moment_z[i]+=mom[2];
}
