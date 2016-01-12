#ifndef UNCERTAINTY_H
#define UNCERTAINTY_H
#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace Eigen;

/*
 Données relatives au Laser
*/
struct Laser{

    //Coordonnées sphériques
    double m_R;
    double m_theta;
    double m_phi;

    //Dérivation des coordonnées sphériques
    Matrix3d m_NablaP;

    //Centre du miroir laser
    Vector3d m_X_o_lidar;

};
/*
 Données relatives aux erreurs Laser (6 sources d'incertitudes)
*/
struct Error_laser{

    // Incertitude les coordonnées r, theta et phi
    double m_delta_R; // metres
    double m_delta_theta;
    double m_delta_phi;

    // Incertitude sur la position du centre du miroir Laser
    double m_delta_beam_origin_x;
    double m_delta_beam_origin_y;
    double m_delta_beam_origin_z;

};

/*
 Données relatives à la centrale inertielle
*/
struct Applanix{

    double m_roll; //rad
    double m_pitch; //rad
    double m_heading; //rad

};
/*
 Données relatives aux erreurs de la centrale inertielle (3 sources d'incertitudes)
*/
struct Error_applanix{

    double m_delta_roll; //  minutes
    double m_delta_pitch; // minutes
    double m_delta_heading; // minutes

};

/*
 Données relatives à la calibration (entre la centrale inertielle et le laser)
*/
struct Calibration{

    //rotation et translation entre laser et centrale inertielle
    Matrix3d m_rotation_calibration;
    Vector3d m_translation_calibration;

};
/*
 Données relatives aux erreurs de calibration (entre la centrale inertielle et le laser).
 (6 sources d'incertitudes)
*/
struct Error_calibration{

    double m_delta_roll_calibration;
    double m_delta_pitch_calibration;
    double m_delta_heading_calibration;

    // translation entre applanix et laser

    double m_Error_translation_x;
    double m_Error_translation_y;
    double m_Error_translation_z;

};
/*
 Données relatives au repère Monde
*/
struct World{

    //Repère Terrain
    double m_longitude;
    double m_latitude;
    double m_altitude;

    //Lambert 93
    double m_conv_meridien;

};
/*
 Données relatives aux erreurs dans le repère Monde (3 sources d'incertitudes)
*/
struct Error_world{

    // Incertitudes sur la longitude, la latitude, et l'altitude
    double m_delta_longitude;
    double m_delta_latitude;
    double m_delta_altitude;


};
/*
 Coordonnées du point dans le repère laser et Applanix
*/
struct Point{
    // Point 3D dans le référentiel Lidar
    Vector3d m_X_lidar;
    // Point 3D dans le référentiel Applanix
    Vector3d m_X_app;


};

/*
La class Uncertainty a pour objectif de calculer la matrice de variance-covariance d'un point du nuage laser.
Pour cela, nous prenons en compte 18 sources d'incertitudes.
*/
class Uncertainty
{
private :

    static double deg_to_rad; // Convertir degrées en radians
    static double arc_minutes_to_degree; // Convertir minutes en degrées
    static Matrix3d Id; // Matrice identité
    static Matrix3d M_NEB_ENH; // Matrice rotation permettant de passer du repere NEB vers ENH

    Laser laser;
    Error_laser error_laser;

    Applanix applanix;
    Error_applanix error_applanix;

    Calibration calibration;
    Error_calibration error_calibration;

    World world;
    Error_world error_world;

    Point point;

public:

    // Première initialisation (attributs)
    void setParam(double longitude,double latitude,double altitude,double roll,double pitch,double heading,double delta_roll,
                                    double delta_pitch,double delta_heading,double conv,double northPositionRMSError,double eastPositionRMSError,double downPositionRMSError);
    // Deuxième initialisation (attributs)
    void setParam(double range,double theta,double phi,const Matrix3d& rotation,const Vector3d& translation,const Vector3d& X_o_Lidar,const Vector3d& X_lidar,const Vector3d& Xapp);

    // Matrice variance-covariance dans le repère terrain
    Matrix3d covariance_terrain();


    double getConvMeridien();
    double getLongitude();
    double getLatitude();
    double getAltitude();
    double getRoll();
    double getPitch();
    double getHeading();



private:


//----------------------------------------------------Méthodes privées de la classe Uncertainty--------------------------------------//



/*Rx: rotation about X-axis, roll
     |1  0   0|
     |0 Cx -Sx|
     |0 Sx  Cx|
*/
inline Matrix3d matrice_Rotation_roll(double roll){

    Matrix3d matrice_rotation;
    double cos_roll=cos(roll),sin_roll=sin(roll);
    matrice_rotation<<1,0,0,
                      0,cos_roll,-sin_roll,
                      0,sin_roll,cos_roll;

    return matrice_rotation;
}


/* Ry: rotation about Y-axis, pitch
          Ry
     | Cy  0 Sy|
     |  0  1  0|
     |-Sy  0 Cy|
*/
inline Matrix3d matrice_Rotation_pitch(double pitch){

    Matrix3d matrice_rotation;
    double cos_pitch=cos(pitch),sin_pitch=sin(pitch);
    matrice_rotation<<cos_pitch,0,sin_pitch,
                      0,1,0,
                      -sin_pitch,0,cos_pitch;

    return matrice_rotation;
}

/* Rz: rotation about Z-axis, Yaw
            Rz
      |Cz -Sz 0|
      |Sz  Cz 0|
      | 0   0 1|
*/
inline Matrix3d matrice_Rotation_heading(double heading){

        Matrix3d matrice_rotation;
        double cos_heading=cos(heading),sin_heading=sin(heading);
        matrice_rotation<<cos_heading,-sin_heading,0,
                          sin_heading,cos_heading,0,
                          0,0,1;

        return matrice_rotation;
}

// ∂Rx/∂roll
inline Matrix3d matrice_Rotation_Droll(double roll){

        Matrix3d matrice_rotation;
        double sin_roll=sin(roll),cos_roll=cos(roll);
        matrice_rotation<<0,0,0,
                          0,-sin_roll,-cos_roll,
                          0,cos_roll,-sin_roll;

        return matrice_rotation;
}

// ∂Ry/∂pitch
inline Matrix3d matrice_Rotation_Dpitch(double pitch){

        Matrix3d matrice_rotation;
        double sin_pitch=sin(pitch),cos_pitch=cos(pitch);

        matrice_rotation<<-sin_pitch,0,cos_pitch,
                          0,0,0,
                          -cos_pitch,0,-sin_pitch;

        return matrice_rotation;
}

// ∂Ry/∂heading
inline Matrix3d matrice_Rotation_Dheading(double heading){

        Matrix3d matrice_rotation;
        double sin_heading=sin(heading),cos_heading=cos(heading);

        matrice_rotation<<-sin_heading,-cos_heading,0,
                          cos_heading,-sin_heading,0,
                          0,0,0;
        return matrice_rotation;

}


/* R rotation about X,Y,Z-axis, roll, pitch and heading
    R=Rz*Ry*Rx =| CzCy        -SzCx+CzSySx         SzSx+CzSyCx  |
                | SzCy         CzCx+SzSySx        -CzSx+SzSyCx  |
                |-Sy           CySx                       CyCx  |
*/
inline Matrix3d matrice_Rotation(double roll,double pitch,double heading){

        Matrix3d matrice_rotation=matrice_Rotation_heading(heading)*matrice_Rotation_pitch(pitch)*matrice_Rotation_roll(roll);

        return matrice_rotation;

}


/*
 Matrice permettant de passer dans le repère Lambert 93
 */
inline Matrix3d Matrice_CONV(double conv){

        Matrix3d mat;
        double cos_conv=cos(conv),sin_conv=sin(conv);

        mat << cos_conv,-sin_conv,0,
               sin_conv,cos_conv,0,
               0,0,1;

        return mat;

}



/*
Determine les 3 angles à partir de la matrice de rotation (entre laser et applanix).
*/

inline Vector3d decompose_rotation(const Matrix3d &R){

       Vector3d decompose_rotation;

       decompose_rotation <<asin(-R(1,2)), atan2(R(0,2), R(2,2)),atan2(R(1,0),R(1,1)); // Formule pour R=Ry.Rz.Rx pour theta_y!=0 et theta_y!= pi

       return decompose_rotation;

}


public:


// Retourne les coordonnées du point du nuage dans le repère terrain

inline Vector3d X_ground(){

    Vector3d angle=decompose_rotation(calibration.m_rotation_calibration);
    Vector3d tran;
    tran<<  world.m_longitude,world.m_latitude, world.m_altitude;

    return Matrice_CONV(world.m_conv_meridien)*M_NEB_ENH*matrice_Rotation(applanix.m_roll,applanix.m_pitch,applanix.m_heading)*
            ( matrice_Rotation_pitch(angle(1))*matrice_Rotation_roll(angle(0))*matrice_Rotation_heading(angle(2))*point.m_X_lidar+calibration.m_translation_calibration)+tran;



}


};

#endif // UNCERTAINTY_H
