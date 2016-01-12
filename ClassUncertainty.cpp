#include "uncertainty.h"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <iomanip>

using namespace Eigen;

namespace {
    //Matrice rotation repère NEB -> ENH
    Matrix3d buildNEB_ENH() {
        Matrix3d NEB_ENH;
        NEB_ENH<<0.,1.,0.,
                1.,0.,0.,
                0.,0.,-1.;
        return NEB_ENH;
    }
}

double Uncertainty::deg_to_rad=3.1415/180.;
double Uncertainty::arc_minutes_to_degree=1./60.;

Matrix3d Uncertainty::Id = Matrix3d::Identity();
Matrix3d Uncertainty::M_NEB_ENH=buildNEB_ENH();


/*Initialisation attributs objet*/
void Uncertainty::setParam(double longitude,double latitude,double altitude,double roll,double pitch,double heading,double delta_roll,
                                 double delta_pitch,double delta_heading,double conv,double northPositionRMSError,double eastPositionRMSError,double downPositionRMSError){
    // Paramètres incertitudes

    world={longitude,latitude,altitude,conv};
    error_world={northPositionRMSError,eastPositionRMSError,downPositionRMSError};


    applanix={roll,pitch,heading};
    error_applanix={delta_roll*arc_minutes_to_degree*deg_to_rad,delta_pitch*arc_minutes_to_degree*deg_to_rad,delta_heading*arc_minutes_to_degree*deg_to_rad};


}


/*Initialisation attributs objet*/
void Uncertainty::setParam(double range,double theta,double phi,const Matrix3d& rotation,const Vector3d& translation,const Vector3d& X_o_Lidar,const Vector3d& X_lidar,const Vector3d& Xapp){

    Matrix3d nablaP;
    double cos_phi=cos(phi),sin_phi=sin(phi),cos_theta=cos(theta),sin_theta=sin(theta);

    nablaP<< -cos_theta*sin_phi,range*sin_theta*sin_phi,-range*cos_theta*cos_phi,
             -sin_theta*sin_phi,-range*cos_theta*sin_phi,-range*sin_theta*cos_phi,
             -cos_phi,0,range*sin_phi;

    // Paramètres incertitudes
    laser={range,theta,phi,nablaP,X_o_Lidar};
    error_laser={0.005,(0.001*deg_to_rad),(0.001*deg_to_rad),0.0001,0.0001,0.0001};

    calibration={rotation,translation};
    error_calibration={0.1*deg_to_rad,0.1*deg_to_rad,0.1*deg_to_rad,0.001,0.001,0.001};

    point={X_lidar,Xapp};


}


/*

 Retourne la matrice de variance-covariance d'un point du nuage laser dans le repère terrain

 Matrice de variance-covariance calculé à partir de la formule de géo-référencement direct suivante :

                                                                                |beam_origin_x + r cos(θ) sin(φ) |
  X_ground= R_conv * R_NEB_NEH * R_Heading * R_Pitch * R_Roll * [R_calibration *|beam_origin_y + r sin(θ) sin(φ) | + T_calibration ] + T_ground
                                                                                |beam_origin_z + r cos(φ)        |

                                                                                |-------------X_lidar------------|
                                                                 |--------------------------Xapp-----------------------------------|
*/

Matrix3d Uncertainty::covariance_terrain(){

    Matrix3d M_CONV=Matrice_CONV(world.m_conv_meridien);
    Matrix3d M_Rotation_HPR=matrice_Rotation(applanix.m_roll,applanix.m_pitch,applanix.m_heading);
    Matrix3d M_Rotation_R=matrice_Rotation_roll(applanix.m_roll);
    Matrix3d M_Rotation_P=matrice_Rotation_pitch(applanix.m_pitch);
    Matrix3d M_Rotation_H=matrice_Rotation_heading(applanix.m_heading);

    //Coordonnées du point dans le repere scanner laser
    const Vector3d& X_lidar=point.m_X_lidar;

    //Coordonnées du point dans le repère Applanix
    const Vector3d& Xapp=point.m_X_app;


    Vector3d angle=decompose_rotation(calibration.m_rotation_calibration);

/*
      Matrice B = ∂F/ ∂l Dérivée par rapport aux observations.
      18 Observations : r,theta,phi,
                    beam_origin_x,beam_origin_y,beam_origin_z,
                    translation_calibration_x,translation_calibration_y,translation_calibration_z,
                    roll_calibration,pitch_calibration,heading_calibration,
                    roll,pitch,heading,
                    translation_applanix_terrain_x,translation_applanix_terrain_y,translation_applanix_terrain_z
*/

    Matrix3d M_terrain=M_CONV*M_NEB_ENH;
    Matrix3d M_terrain_HPR_calibration=M_terrain*M_Rotation_HPR*calibration.m_rotation_calibration;
    Vector3d M_Rotation_R_Xapp=M_Rotation_R*Xapp;


    MatrixXd B(3,18),vec(18,1);
    B<<      M_terrain_HPR_calibration*laser.m_NablaP,
            -M_terrain_HPR_calibration,
            -M_terrain*M_Rotation_HPR,
            -M_terrain*M_Rotation_HPR*matrice_Rotation_pitch(angle(1))*matrice_Rotation_Droll(angle(0))*matrice_Rotation_heading(angle(2))*X_lidar,
            -M_terrain*M_Rotation_HPR*matrice_Rotation_Dpitch(angle(1))*matrice_Rotation_roll(angle(0))*matrice_Rotation_heading(angle(2))*X_lidar,
            -M_terrain*M_Rotation_HPR*matrice_Rotation_pitch(angle(1))*matrice_Rotation_roll(angle(0))*matrice_Rotation_Dheading(angle(2))*X_lidar,
            -M_terrain*M_Rotation_H*M_Rotation_P*matrice_Rotation_Droll(applanix.m_roll)*Xapp,
            -M_terrain*M_Rotation_H*matrice_Rotation_Dpitch(applanix.m_pitch)*M_Rotation_R_Xapp,
            -M_terrain*matrice_Rotation_Dheading(applanix.m_heading)*M_Rotation_P*M_Rotation_R_Xapp,
            -M_terrain*Id;

    //Vecteur contenant les 18 sources d'incertitudes
    vec<<error_laser.m_delta_R*error_laser.m_delta_R,
         error_laser.m_delta_theta*error_laser.m_delta_theta,
         error_laser.m_delta_phi*error_laser.m_delta_phi,

         error_laser.m_delta_beam_origin_x*error_laser.m_delta_beam_origin_x,
         error_laser.m_delta_beam_origin_y*error_laser.m_delta_beam_origin_y,
         error_laser.m_delta_beam_origin_z*error_laser.m_delta_beam_origin_z,

         error_calibration.m_Error_translation_x*error_calibration.m_Error_translation_x,
         error_calibration.m_Error_translation_y*error_calibration.m_Error_translation_y,
         error_calibration.m_Error_translation_z*error_calibration.m_Error_translation_z,

         error_calibration.m_delta_roll_calibration*error_calibration.m_delta_roll_calibration,
         error_calibration.m_delta_pitch_calibration*error_calibration.m_delta_pitch_calibration,
         error_calibration.m_delta_heading_calibration*error_calibration.m_delta_heading_calibration,

         error_applanix.m_delta_roll*error_applanix.m_delta_roll,
         error_applanix.m_delta_pitch*error_applanix.m_delta_pitch,
         error_applanix.m_delta_heading*error_applanix.m_delta_heading,

         error_world.m_delta_longitude*error_world.m_delta_longitude,
         error_world.m_delta_latitude*error_world.m_delta_latitude,
         error_world.m_delta_altitude*error_world.m_delta_altitude;

      MatrixXd B_dot_vec(3,18); // Produit matriciel de B et vec

      for (int i = 0; i< B.rows(); i++){
              for (int j = 0; j< vec.size(); j++){ B_dot_vec(i,j)=B(i,j)*vec(j);}
      }

      return B_dot_vec*B.transpose();

}



double Uncertainty::getConvMeridien(){

    return world.m_conv_meridien;

}

double Uncertainty::getLongitude(){

    return world.m_longitude;

}

double Uncertainty::getLatitude(){

    return world.m_latitude;

}
double Uncertainty::getAltitude(){

    return world.m_altitude;

}

double Uncertainty::getRoll(){

    return applanix.m_roll;
}

double Uncertainty::getPitch(){

    return applanix.m_pitch;
}

double Uncertainty::getHeading(){

    return applanix.m_heading;
}

