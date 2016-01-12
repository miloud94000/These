
#include "XMlsUnc.h"
//#include "libXBaseXML/XArchiXMLException.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

//namespace {
    //Matrice rotation repère NEB -> ENH
    Matrix3d buildNEB_ENH() {
        Matrix3d NEB_ENH;
        NEB_ENH<<0.,1.,0.,
                1.,0.,0.,
                0.,0.,-1.;
        return NEB_ENH;
    }
//}

double deg_to_rad=3.1415/180.;
double arc_minutes_to_degree=1./60.;

Matrix3d Id = Matrix3d::Identity();
Matrix3d M_NEB_ENH=buildNEB_ENH();

/*

 Retourne la matrice de variance-covariance d'un point du nuage laser

 Matrice de variance-covariance calculé à partir de la formule de géo-référencement direct suivante :

                                                                                |beam_origin_x + r cos(θ) cos(φ) |
  X_ground= R_conv * R_NEB_NEH * R_Heading * R_Pitch * R_Roll * [R_calibration *|beam_origin_y + r sin(θ) cos(φ) | + T_calibration ] + T_ground
                                                                                |beam_origin_z + r sin(φ)        |

                                                                                |-------------X_Laser------------|
                                                                 |--------------------------X_applanix-----------------------------|
Warning : beam_origin n'est pas visible par l'utilisateur
*/

Matrix3d XMlsUnc::Uncertainty(XEchoIndex i_echo, XVariance variance)
{
        float range=Range(i_echo);
        float theta=Theta(i_echo);
        float phi= Phi(i_echo); // Attention l'angle phi ne corespond pas à la convention des coordonnées sphériques

        double cos_phi=cos(phi),sin_phi=sin(phi),cos_theta=cos(theta),sin_theta=sin(theta);

        // Coordonnées du point dans le repère Laser
        Vector3d X_Laser;
        X_Laser << range*cos_theta*cos_phi ,range*sin_theta*cos_phi, range*sin_phi;

        // Matrice de rotation Laser-> INS (Calibration)
        XMat3D  R=m_sensor_georef.Rotation();

        Matrix3d R_calibration;
        R_calibration << R.A.X,R.A.Y,R.A.Z,
                         R.B.X,R.B.Y,R.B.Z,
                         R.C.X,R.C.Y,R.C.Z;

       // Vecteur translation Laser-> INS (Calibration)
       Vector3d T_calibration;
       T_calibration<< m_sensor_georef.Translation().X , m_sensor_georef.Translation().Y,m_sensor_georef.Translation().Z;

       //Coordonnées du point dans le repère de la centrale inertielle
       Vector3d X_applanix;
       X_applanix=R_calibration*X_Laser+T_calibration;

       // Paramètres de l'INS
       SbetEvent paramINS=Sbet(i_echo);

       // Convergence des méridiens
       double conv_meridien=-m_trajecto.SbetSeries().GetConvMeridien(paramINS);

       // Incertitudes de la centrale inertielle
       AccuracyEvent incertitudeINS= Accuracy(i_echo); // incertitudes INS

        double sigmaRoll=incertitudeINS.m_RollRMSError*arc_minutes_to_degree*deg_to_rad;
        double sigmaPitch=incertitudeINS.m_PitchRMSError*arc_minutes_to_degree*deg_to_rad;
        double sigmaHeading=incertitudeINS.m_headingRMSError*arc_minutes_to_degree*deg_to_rad;

        // Translation INS-> terrain (Lambert 93)
        Vector3d translation_applanix_terrain;
        XArchiGeoref ins =  Ins(i_echo, false);

        translation_applanix_terrain<< ins.Translation().X,ins.Translation().Y,ins.Translation().Z;

        // Matrice 3*3 [∂X_laser/∂r ∂X_laser/∂theta ∂X_laser/∂phi]
        Matrix3d nablaP;
        nablaP<< -cos_theta*cos_phi, range*sin_theta*cos_phi,range*cos_theta*sin_phi,
                 - sin_theta*cos_phi,-range*cos_theta*cos_phi,range*sin_theta*sin_phi,
                 - sin_phi,0,-range*cos_phi;

        Matrix3d M_CONV=Matrice_CONV(conv_meridien); // Matrice rotation terrain (ENH)-> terrain (Lambert 93)
        Matrix3d M_Rotation_HPR=matrice_Rotation(paramINS.m_roll ,paramINS.m_pitch,paramINS.m_plateformHeading); // Matrice rotation INS->terrain ENB
        Matrix3d M_Rotation_R=matrice_Rotation_roll(paramINS.m_roll); // Rotation autour de l'axe X (INS->terrain (ENB))
        Matrix3d M_Rotation_P=matrice_Rotation_pitch(paramINS.m_pitch); // Rotation autour de l'axe y (INS->terrain (ENB))
        Matrix3d M_Rotation_H=matrice_Rotation_heading(paramINS.m_plateformHeading); //Rotation autour de l'axe z (INS->terrain (ENB))

        Vector3d angle=decompose_rotation(R_calibration);



        Matrix3d M_terrain=M_CONV*M_NEB_ENH;
        Matrix3d M_terrain_HPR_calibration=M_terrain*M_Rotation_HPR*R_calibration;
        Vector3d M_Rotation_R_Xapp=M_Rotation_R*X_applanix;


        //     Matrice B = ∂F/ ∂l Dérivée par rapport aux observations.
        //      18 Observations : r,theta,phi, (attention à l'angle phi ne corespond pas au phi des coordonnées sphériques)
        //                    beam_origin_x,beam_origin_y,beam_origin_z, (Position du miroir Laser)
        //                    translation_calibration_x,translation_calibration_y,translation_calibration_z,
        //                    roll_calibration,pitch_calibration,heading_calibration, (rotation calibration)
        //                    roll,pitch,heading, (centrale inertielle)
        //                    translation_applanix_terrain_x,translation_applanix_terrain_y,translation_applanix_terrain_z


        MatrixXd B(3,18);
        B<<      M_terrain_HPR_calibration*nablaP,
                -M_terrain_HPR_calibration,
                -M_terrain*M_Rotation_HPR,
                -M_terrain*M_Rotation_HPR*matrice_Rotation_pitch(angle(1))*matrice_Rotation_Droll(angle(0))*matrice_Rotation_heading(angle(2))*X_Laser,
                -M_terrain*M_Rotation_HPR*matrice_Rotation_Dpitch(angle(1))*matrice_Rotation_roll(angle(0))*matrice_Rotation_heading(angle(2))*X_Laser,
                -M_terrain*M_Rotation_HPR*matrice_Rotation_pitch(angle(1))*matrice_Rotation_roll(angle(0))*matrice_Rotation_Dheading(angle(2))*X_Laser,
                -M_terrain*M_Rotation_H*M_Rotation_P*matrice_Rotation_Droll(paramINS.m_roll)*X_applanix,
                -M_terrain*M_Rotation_H*matrice_Rotation_Dpitch(paramINS.m_pitch)*M_Rotation_R_Xapp,
                -M_terrain*matrice_Rotation_Dheading(paramINS.m_plateformHeading)*M_Rotation_P*M_Rotation_R_Xapp,
                -M_terrain*Id;


       //Vecteur contenant les 18 sources d'incertitudes
        MatrixXd vec(18,1);
        vec<<variance.sigmaR*variance.sigmaR,
              variance.sigmaTheta*variance.sigmaTheta,
              variance.sigmaPhi*variance.sigmaPhi,

              variance.sigmaOx*variance.sigmaOx,
              variance.sigmaOy*variance.sigmaOy,
              variance.sigmaOz*variance.sigmaOz,

              variance.sigmaTranslationCalibrationx*variance.sigmaTranslationCalibrationx,
              variance.sigmaTranslationCalibrationy*variance.sigmaTranslationCalibrationy,
              variance.sigmaTranslationCalibrationz*variance.sigmaTranslationCalibrationz,

              variance.sigmaRollCalibration*variance.sigmaRollCalibration,
              variance.sigmaPitchCalibration*variance.sigmaPitchCalibration,
              variance.sigmaHeadingCalibration*variance.sigmaHeadingCalibration,

              sigmaRoll*sigmaRoll,
              sigmaPitch*sigmaPitch,
              sigmaHeading*sigmaHeading,

              incertitudeINS.m_northPositionRMSError*incertitudeINS.m_northPositionRMSError,
              incertitudeINS.m_eastPositionRMSError*incertitudeINS.m_eastPositionRMSError,
              incertitudeINS.m_downPositionRMSError*incertitudeINS.m_downPositionRMSError;

          MatrixXd B_dot_vec(3,18); // Produit matriciel de B et vec

          for (int i = 0; i< B.rows(); i++){
                  for (int j = 0; j< vec.size(); j++){ B_dot_vec(i,j)=B(i,j)*vec(j);}
          }

          return B_dot_vec*B.transpose();
}
