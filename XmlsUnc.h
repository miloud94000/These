#pragma once

/// Mobile laser scan class with uncertainties

#include "XMls.h"
#include <eigen3/Eigen/Dense>
//using namespace Eigen;

/// Uncertainties on all variables leading to computation of a point geometry
struct XVariance
{
    double sigmaPhi, sigmaTheta, sigmaR;
    double sigmaOx, sigmaOy, sigmaOz;
    double sigmaTranslationCalibrationx,sigmaTranslationCalibrationy,sigmaTranslationCalibrationz;
    double sigmaRollCalibration,sigmaPitchCalibration,sigmaHeadingCalibration;

    XVariance():sigmaPhi(0.001*3.1415/180.), sigmaTheta(0.001*3.1415/180.), sigmaR(0.005),
        sigmaOx(0.001), sigmaOy(0.001), sigmaOz(0.001),
        sigmaTranslationCalibrationx(0.001),sigmaTranslationCalibrationy(0.001),sigmaTranslationCalibrationz(0.001),
        sigmaRollCalibration(0.1*3.1415/180.),sigmaPitchCalibration(0.1*3.1415/180.),sigmaHeadingCalibration(0.1*3.1415/180.){}
    // etc
};

class XMlsUnc:public XMls
{
public:
    XMlsUnc(std::string ept_folder, std::string laser_calib,
         std::string sbet_folder, std::string sbet_mission):
        XMls(ept_folder,laser_calib,sbet_folder,sbet_mission){}

    /// Variace covariance matrix on position of Pworld(XEchoIndex i_echo)
    Eigen::Matrix3d Uncertainty(XEchoIndex i_echo, XVariance variance=XVariance());


private:


    /*Rx: rotation about X-axis, roll
     |1  0   0|
     |0 Cx -Sx|
     |0 Sx  Cx|
*/
    inline Eigen::Matrix3d matrice_Rotation_roll(double roll){

        Eigen::Matrix3d matrice_rotation;
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
    inline Eigen::Matrix3d matrice_Rotation_pitch(double pitch){

        Eigen::Matrix3d matrice_rotation;
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
    inline Eigen::Matrix3d matrice_Rotation_heading(double heading){

        Eigen::Matrix3d matrice_rotation;
        double cos_heading=cos(heading),sin_heading=sin(heading);
        matrice_rotation<<cos_heading,-sin_heading,0,
                sin_heading,cos_heading,0,
                0,0,1;

        return matrice_rotation;
    }

    // ∂Rx/∂roll
    inline Eigen::Matrix3d matrice_Rotation_Droll(double roll){

        Eigen::Matrix3d matrice_rotation;
        double sin_roll=sin(roll),cos_roll=cos(roll);
        matrice_rotation<<0,0,0,
                0,-sin_roll,-cos_roll,
                0,cos_roll,-sin_roll;

        return matrice_rotation;
    }

    // ∂Ry/∂pitch
    inline Eigen::Matrix3d matrice_Rotation_Dpitch(double pitch){

        Eigen::Matrix3d matrice_rotation;
        double sin_pitch=sin(pitch),cos_pitch=cos(pitch);

        matrice_rotation<<-sin_pitch,0,cos_pitch,
                0,0,0,
                -cos_pitch,0,-sin_pitch;

        return matrice_rotation;
    }

    // ∂Ry/∂heading
    inline Eigen::Matrix3d matrice_Rotation_Dheading(double heading){

        Eigen::Matrix3d matrice_rotation;
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
    inline Eigen::Matrix3d matrice_Rotation(double roll,double pitch,double heading){

        Eigen::Matrix3d matrice_rotation=matrice_Rotation_heading(heading)*matrice_Rotation_pitch(pitch)*matrice_Rotation_roll(roll);

        return matrice_rotation;

    }


    /*
 Matrice permettant de passer dans le repère Lambert 93
 */
    inline Eigen::Matrix3d Matrice_CONV(double conv){

        Eigen::Matrix3d mat;
        double cos_conv=cos(conv),sin_conv=sin(conv);

        mat << cos_conv,-sin_conv,0,
                sin_conv,cos_conv,0,
                0,0,1;

        return mat;

    }



    /*
Determine les 3 angles à partir de la matrice de rotation R=Ry(angle0).Rx(angle1).Rz(angle2) (entre laser et applanix).
*/

    inline Eigen::Vector3d decompose_rotation(const Eigen::Matrix3d &R){

        Eigen::Vector3d decompose_rotation;

        decompose_rotation <<asin(-R(1,2)), atan2(R(0,2), R(2,2)),atan2(R(1,0),R(1,1)); // Formule pour R=Ry.Rx.Rz pour theta_y!=0 et theta_y!= pi

        return decompose_rotation;

    }


};
