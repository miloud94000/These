#ifndef __COVARIANCEMATRIX_H__
#define __COVARIANCEMATRIX_H__

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <functional>
#include <numeric>
#include "ReaderFilePosPac/trajecto_reader.h"

#include <Eigen/Sparse>
#include <lib3ds/file.h>
#include <lib3ds/mesh.h>
#include <lib3ds/node.h>
#include <lib3ds/material.h>
#include "LgPoint3.hpp"
#include "ClassKdTree.h"
#include "ClassTriangle.h"
#include "accelerators/kdtreeaccel.h"
//#include "core/pbrt.h"
//#include "core/imageio.cpp"

class Covariance
{

/*!
 *   Matrice de variance-covariance de la position du véhicule (ou trajectoire) du véhicule du temps T0 à Tn (à 1Hz)
 *
 *
 !*/
public:


     Covariance(Trajecto & trajecto,const double & amplification_x=1.,const double & amplification_y=1.,const double & amplification_z=1.); // à 1Hz


     Eigen::MatrixXd const& getCovarianceX(){return m_covarianceX;}    // Matrice variance-covariance trajectoire 1D x  à 1Hz
     Eigen::MatrixXd const& getCovarianceY(){return m_covarianceY;}    // Matrice variance-covariance trajectoire 1D y  à 1Hz
     Eigen::MatrixXd const& getCovarianceZ(){return m_covarianceZ;}    // Matrice variance-covariance trajectoire 1D z  à 1Hz
     Eigen::MatrixXd const& getCovariance(){return m_covariance;}      // Matrice variance-covariance trajectoire 3D xyz  à 1Hz

     Eigen::MatrixXd getCovarianceXY(){return m_covarianceXY;}  // Matrice variance-covariance trajectoire 2D xy  à 1Hz


private :

    TEventSeries<AccuracyEvent> & m_AccEvent; // Incertitude sur la trajectoire
    int m_NbreEvenement;                      // Nombre d'événement sur la trajectoire du véhicule
    Trajecto m_trajecto;

    Eigen::MatrixXd m_covarianceX; // Covariance 1D X
    /*!      t0x t1x ... tnx
     *  t0x [              ]
     *      [              ]
     *      [              ]
     *  tnx [              ]
     * !*/

    Eigen::MatrixXd m_covarianceY;   // Covariance 1D Y
    /*!      t0y t1y ... tny
     *  t0y [              ]
     *      [              ]
     *      [              ]
     *  tny [              ]
     * !*/

    Eigen::MatrixXd m_covarianceZ;   // Covariance 1D Z
    /*!      t0z t1z ... tnz
     *  t0z [              ]
     *      [              ]
     *      [              ]
     *  tnz [              ]
     * !*/


   Eigen::MatrixXd m_covarianceXY;   // Covariance 2D Warning : à 1 Hz (Necessaire pour le recalage)
    /*!      t0x t0y t1x t1y... tnx tny
     *  t0x [                          ]
     *   .  [                          ]
     *  tnx [                          ]
     *  tny [                          ]
     * !*/


    Eigen::MatrixXd m_covariance;    // Covariance 3D Warning : à 1Hz
    /*!      t0x t0y t0z t1x t1y t1z... tnx tny tnz
     *  t0x [                                      ]
     *   .  [                                      ]
     *    . [                                      ]
     *  tnz [                                      ]
     * !*/





    inline double facteur_normalisation(const double & sigma_p_i,const double & sigma_p_i_1,const double  & sigma_v_i,const double & dt){

         return sqrt((sigma_p_i_1*sigma_p_i_1)/ ((dt*dt*sigma_v_i*sigma_v_i)+(sigma_p_i*sigma_p_i)) );
    }


    inline Lg::Point3d Propagation_uncertainty_from_NEB_to_LAM93_Position(const SbetEvent&  sbetEvent, const AccuracyEvent & accEvent){


        Eigen::Matrix3d Position_uncertainty_vehicule;
        Position_uncertainty_vehicule << accEvent.m_northPositionRMSError*accEvent.m_northPositionRMSError,0,0,
                                        0,accEvent.m_eastPositionRMSError*accEvent.m_eastPositionRMSError,0,
                                        0,0,accEvent.m_downPositionRMSError*accEvent.m_downPositionRMSError; // Matrice covariance de la position du vehicule (LAMB93)

       Eigen::Matrix3d NEB_ENH;
       NEB_ENH << 0.,1.,0.,
                  1.,0.,0.,
                  0.,0.,-1.;

      const double & conv =-m_trajecto.SbetSeries().GetConvMeridien(sbetEvent);
      double cos_conv=cos(conv),sin_conv=sin(conv);

       Eigen::Matrix3d Rconv;
       Rconv << cos_conv,-sin_conv,0,
                sin_conv,cos_conv,0,
                0,0,1;

       Eigen::Matrix3d A;

       A=Rconv*NEB_ENH;
       Position_uncertainty_vehicule=A*Position_uncertainty_vehicule*A.transpose();

       Lg::Point3d Pos_uncertainty_lamb93(sqrt(Position_uncertainty_vehicule(0,0)),sqrt(Position_uncertainty_vehicule(1,1)),sqrt(Position_uncertainty_vehicule(2,2)));


       return Pos_uncertainty_lamb93; // Incertitude dans le repère Lamb93


    }


    inline Lg::Point3d Propagation_uncertainty_from_NEB_to_LAM93_Vitesse(const SbetEvent&  sbetEvent, const AccuracyEvent & accEvent){


        Eigen::Matrix3d Position_uncertainty_vehicule;
        Position_uncertainty_vehicule << accEvent.m_northVelocityRMSError*accEvent.m_northVelocityRMSError,0,0,
                                        0,accEvent.m_eastVelocityRMSError*accEvent.m_eastVelocityRMSError,0,
                                        0,0,accEvent.m_downVelocityRMSError*accEvent.m_downVelocityRMSError; // Matrice covariance de la position du vehicule (LAMB93)

       Eigen::Matrix3d NEB_ENH;
       NEB_ENH<<0.,1.,0.,
                1.,0.,0.,
                0.,0.,-1.;

      const double & conv =-m_trajecto.SbetSeries().GetConvMeridien(sbetEvent);
      double cos_conv=cos(conv),sin_conv=sin(conv);

       Eigen::Matrix3d Rconv;
       Rconv << cos_conv,-sin_conv,0,
                sin_conv,cos_conv,0,
                0,0,1;

       Eigen::Matrix3d A;

       A=Rconv*NEB_ENH;
       Position_uncertainty_vehicule=A*Position_uncertainty_vehicule*A.transpose();

       Lg::Point3d Pos_uncertainty_lamb93Vitesse(sqrt(Position_uncertainty_vehicule(0,0)),sqrt(Position_uncertainty_vehicule(1,1)),sqrt(Position_uncertainty_vehicule(2,2)));

       return Pos_uncertainty_lamb93Vitesse; // Incertitude dans le repère Lamb93


    }


};



/*
 *
 *
 *   Pour retourner un vecteur aléatoire suivant une gaussienne . EN entrée une covariance et une moyenne (vecteur).
 *
 *
 */
namespace Eigen {
  namespace internal {
    template<typename Scalar>
      struct scalar_normal_dist_op
      {
    static std::mt19937 rng;                        // The uniform pseudo-random algorithm
    mutable std::normal_distribution<Scalar> norm; // gaussian combinator

    EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

    template<typename Index>
    inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
    inline void seed(const uint64_t &s) { rng.seed(s); }
      };

    template<typename Scalar>
      std::mt19937 scalar_normal_dist_op<Scalar>::rng;

    template<typename Scalar>
      struct functor_traits<scalar_normal_dist_op<Scalar> >
      { enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };

  } // end namespace internal

  /**
    Find the eigen-decomposition of the covariance matrix
    and then store it for sampling from a multi-variate normal
  */
  template<typename Scalar>
    class EigenMultivariateNormal
  {
    Matrix<Scalar,Dynamic,Dynamic> _covar;
    Matrix<Scalar,Dynamic,Dynamic> _transform;
    Matrix< Scalar, Dynamic, 1> _mean;
    internal::scalar_normal_dist_op<Scalar> randN; // Gaussian functor
    bool _use_cholesky;
    SelfAdjointEigenSolver<Matrix<Scalar,Dynamic,Dynamic> > _eigenSolver; // drawback: this creates a useless eigenSolver when using Cholesky decomposition, but it yields access to eigenvalues and vectors

  public:
  EigenMultivariateNormal(const Matrix<Scalar,Dynamic,1>& mean,const Matrix<Scalar,Dynamic,Dynamic>& covar,
              const bool use_cholesky=false,const uint64_t &seed=std::mt19937::default_seed)
      :_use_cholesky(use_cholesky)
     {
        randN.seed(seed);
    setMean(mean);
    setCovar(covar);
      }

    void setMean(const Matrix<Scalar,Dynamic,1>& mean) { _mean = mean; }
    void setCovar(const Matrix<Scalar,Dynamic,Dynamic>& covar)
    {
      _covar = covar;

      // Assuming that we'll be using this repeatedly,
      // compute the transformation matrix that will
      // be applied to unit-variance independent normals

      if (_use_cholesky)
    {
      Eigen::LLT<Eigen::Matrix<Scalar,Dynamic,Dynamic> > cholSolver(_covar);
      // We can only use the cholesky decomposition if
      // the covariance matrix is symmetric, pos-definite.
      // But a covariance matrix might be pos-semi-definite.
      // In that case, we'll go to an EigenSolver
      if (cholSolver.info()==Eigen::Success)
        {
          // Use cholesky solver
          _transform = cholSolver.matrixL();
        }
      else
        {
          throw std::runtime_error("Failed computing the Cholesky decomposition. Use solver instead");
        }
    }
      else
    {
      _eigenSolver = SelfAdjointEigenSolver<Matrix<Scalar,Dynamic,Dynamic> >(_covar);
      _transform = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
    }
    }

    /// Draw nn samples from the gaussian and return them
    /// as columns in a Dynamic by nn matrix
    Matrix<Scalar,Dynamic,-1> samples(int nn)
      {
    return (_transform * Matrix<Scalar,Dynamic,-1>::NullaryExpr(_covar.rows(),nn,randN)).colwise() + _mean;
      }
  }; // end class EigenMultivariateNormal
} // end namespace Eigen



#endif


#include "CovarianceMatrix.h"

template <typename T> inline constexpr T carre(T const& v) { return v * v; }

// Calcul la matrice de covariance à 1Hz  tox toy toz ....tnx tny tnz en 3D

// Warning : Covariance dans le repère NEH
Covariance::Covariance(Trajecto &  trajecto,const double & amplification_x,const double & amplification_y,const double & amplification_z) :
          m_AccEvent(trajecto.AccuracySeries()),
          m_NbreEvenement(m_AccEvent.Nevent()),
          m_covarianceX(Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement)),
          m_covarianceY(Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement)),
          m_covarianceZ(Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement)),
          m_covarianceXY(Eigen::MatrixXd::Zero(m_NbreEvenement*2,m_NbreEvenement*2)),
          m_covariance(Eigen::MatrixXd::Zero(m_NbreEvenement*3, m_NbreEvenement*3))
{




    std::vector<double> v_normaliseX;
    std::vector<double> v_normaliseY;
    std::vector<double> v_normaliseZ;
    double para=0.5;


    assert(m_NbreEvenement > 0);
    for(int k=0 ; k<(m_NbreEvenement-1) ; ++k){

        auto const& ek1 = m_AccEvent.Event(k+1);
        auto const& ek = m_AccEvent.Event(k);

        const double dt=ek1.m_time - ek.m_time;

        v_normaliseX.push_back(facteur_normalisation(amplification_y*para/*ek.m_eastPositionRMSError*/,amplification_y*para/*ek1.m_eastPositionRMSError*/,0.007/*ek.m_eastVelocityRMSError*/,dt));
        v_normaliseY.push_back(facteur_normalisation(amplification_x*para/*ek.m_northPositionRMSError*/,amplification_x*para/*ek1.m_northPositionRMSError*/,0.007/*ek.m_northVelocityRMSError*/,dt));
        v_normaliseZ.push_back(facteur_normalisation(amplification_z*ek1.m_downPositionRMSError,amplification_z*ek1.m_downPositionRMSError,ek.m_downVelocityRMSError,dt));

    }

    for(int i=0 ; i<m_NbreEvenement; ++i){

        auto const& event_i = m_AccEvent.Event(i);

        const double & cxi = carre(amplification_y*para/*event_i.m_eastPositionRMSError*/);
        const double & cyi = carre(amplification_x*para/*event_i.m_northPositionRMSError*/);
        const double & czi = carre(amplification_z*event_i.m_downPositionRMSError);

        m_covarianceX(i,i) = cxi; // Covariance 1D x
        m_covarianceY(i,i) = cyi; // Covariance 1D y
        m_covarianceZ(i,i) = czi; // Covariance 1D z

        double normX =1.;
        double normY = 1.;
        double normZ= 1.;

        for(int j=i+1 ; j<m_NbreEvenement; ++j){
              // Construction matrice m_covarianceX
            normX *= v_normaliseX[j-1];
            m_covarianceX(i,j) = cxi *normX;
            m_covarianceX(j,i) = m_covarianceX(i,j);
              // Construction matrice m_covarianceY
            normY *= v_normaliseY[j-1];
            m_covarianceY(i,j) = cyi *normY ;
            m_covarianceY(j,i) = m_covarianceY(i,j);

             // Construction matrice m_covarianceZ
            normZ*=v_normaliseZ[j-1];
            m_covarianceZ(i,j) = czi *normZ  ;
            m_covarianceZ(j,i)=m_covarianceZ(i,j);

        }
   }


    // Construction matrice m_covariance et m_covarianceXY à l'aide des matrice 1D en x ,y et z
     for(int i=0 ; i<m_NbreEvenement; ++i){
           for(int j=i; j<m_NbreEvenement; ++j){

               const double & Cx_ij=m_covarianceX(i,j);
               const double & Cy_ij=m_covarianceY(i,j);
               const double & Cz_ij=m_covarianceZ(i,j);

               const unsigned int & k=i*3;
               const unsigned int & p=j*3;

               const unsigned int & I=i*2;
               const unsigned int & J=j*2;

               // covariance 2D
               m_covarianceXY(I,J)=Cx_ij;
               m_covarianceXY(I+1,J+1)=Cy_ij;


               // covariance 3D
               m_covariance(k,p)=Cx_ij;
               m_covariance(k+1,p+1)=Cy_ij;
               m_covariance(k+2,p+2)=Cz_ij;

               if(i!=j){ // Matrice symétique

                   m_covarianceXY(J,I)=Cx_ij;
                   m_covarianceXY(J+1,I+1)=Cy_ij;


                   m_covariance(p,k)=Cx_ij;
                   m_covariance(p+1,k+1)=Cy_ij;
                   m_covariance(p+2,k+2)=Cz_ij;
               }

         }
    }


}


//// Calcul la matrice de covariance de  la trajectoire à 100Hz (A chaque event de la trajectoire)
//Covariance::Covariance(Trajecto &  trajecto,const double & amplification_x,const double & amplification_y,const double & amplification_z) :
//          m_trajecto(trajecto),
//          m_Event(trajecto.SbetSeries()),
//          m_NbreEvenement(trajecto.SbetSeries().Nevent())
//{


//    clock_t start = clock();
//    std::vector<double> v_normaliseX;
//    std::vector<double> v_normaliseY;
////    std::vector<double> v_normaliseZ;

//    assert(m_NbreEvenement > 0);
//    for(unsigned int k=0 ; k<(m_NbreEvenement-1) ; ++k){

//        const double & time=m_Event.Event(k).m_time;
//        const int & Acc_evt= trajecto.AccuracySeries().PrevIndex(time);
//        AccuracyEvent & ek;
//        trajecto.AccuracySeries().Interpol_event(ek,time);
//        //=Interpol(trajecto.AccuracySeries().Event(Acc_evt), trajecto.AccuracySeries().Event(Acc_evt+1) ,time); // Incertitude au temps t evenement k

//        const double & time1=m_Event.Event(k+1).m_time;
//        const int &  Acc_evt1= trajecto.AccuracySeries().PrevIndex(time1);
//        const AccuracyEvent & ek1=Interpol(trajecto.AccuracySeries().Event(Acc_evt1), trajecto.AccuracySeries().Event(Acc_evt1+1) ,time1) ; // Incertitude au temps t+1   evenement k+1

//        const double & dt=time1 - time;

//         Lg::Point3d Uncertainty_pos_ek=Propagation_uncertainty_from_NEB_to_LAM93_Position(m_Event.Event(k),ek);
//         Lg::Point3d Uncertainty_velocity_ek=Propagation_uncertainty_from_NEB_to_LAM93_Vitesse(m_Event.Event(k),ek);
//         Lg::Point3d Uncertainty_pos_ek1=Propagation_uncertainty_from_NEB_to_LAM93_Position(m_Event.Event(k+1),ek1); // Incertitude de vitesse de l'evenement k


//        v_normaliseX.push_back( facteur_normalisation(amplification_x*Uncertainty_pos_ek.x(),amplification_x*Uncertainty_pos_ek1.x(), Uncertainty_velocity_ek.x(),dt) );
//        v_normaliseY.push_back( facteur_normalisation(amplification_y*Uncertainty_pos_ek.y(),amplification_y*Uncertainty_pos_ek1.y(), Uncertainty_velocity_ek.y(),dt) );
////        v_normaliseZ.push_back( facteur_normalisation(amplification_z*Uncertainty_pos_ek.z(),amplification_z*Uncertainty_pos_ek1.z(), Uncertainty_velocity_ek.z(),dt) );

//    }

//       m_covarianceX=Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement);
//       m_covarianceY=Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement);
////       m_covarianceZ=Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement);

//    for(unsigned int i=0 ; i<m_NbreEvenement; ++i){

//        auto const& event_i = m_Event.Event(i);
//        const double & temps=event_i.m_time;
//        const int &  Accuracy_evt= trajecto.AccuracySeries().PrevIndex(temps);
//        auto const& ek=Interpol(trajecto.AccuracySeries().Event(Accuracy_evt), trajecto.AccuracySeries().Event(Accuracy_evt) ,temps) ; // Récupère l'incertitude au temps t

//        Lg::Point3d Uncertainty_pos_event=Propagation_uncertainty_from_NEB_to_LAM93_Position(m_Event.Event(i),ek);

//        const double & cxi = carre(amplification_x*Uncertainty_pos_event.x());
//        const double & cyi = carre(amplification_y*Uncertainty_pos_event.y());
////        const double czi = carre(amplification_z*Uncertainty_pos_event.z());

//        m_covarianceX(i,i) = cxi;
//        m_covarianceY(i,i) = cyi;
////        m_covarianceZ(i,i) = czi;

//        double normX =1.;
//        double normY = 1.;
////        double normZ= 1.;

//        for(unsigned int j=i+1 ; j<m_NbreEvenement; ++j){
//              // Construction matrice m_covarianceX
//            normX *= v_normaliseX[j-1];
//            m_covarianceX(i,j) = cxi *normX;
//            m_covarianceX(j,i) = m_covarianceX(i,j);
//              // Construction matrice m_covarianceY
//            normY *= v_normaliseY[j-1];
//            m_covarianceY(i,j) = cyi *normY ;
//            m_covarianceY(j,i) = m_covarianceY(i,j);

//// Calcul la matrice de covariance de  la trajectoire à 100Hz (A chaque event de la trajectoire)
//Covariance::Covariance(Trajecto &  trajecto,const double & amplification_x,const double & amplification_y,const double & amplification_z) :
//          m_trajecto(trajecto),
//          m_Event(trajecto.SbetSeries()),
//          m_NbreEvenement(trajecto.SbetSeries().Nevent())
//{


//    clock_t start = clock();
//    std::vector<double> v_normaliseX;
//    std::vector<double> v_normaliseY;
////    std::vector<double> v_normaliseZ;

//    assert(m_NbreEvenement > 0);
//    for(unsigned int k=0 ; k<(m_NbreEvenement-1) ; ++k){

//        const double & time=m_Event.Event(k).m_time;
//        const int & Acc_evt= trajecto.AccuracySeries().PrevIndex(time);
//        AccuracyEvent & ek;
//        trajecto.AccuracySeries().Interpol_event(ek,time);
//        //=Interpol(trajecto.AccuracySeries().Event(Acc_evt), trajecto.AccuracySeries().Event(Acc_evt+1) ,time); // Incertitude au temps t evenement k

//        const double & time1=m_Event.Event(k+1).m_time;
//        const int &  Acc_evt1= trajecto.AccuracySeries().PrevIndex(time1);
//        const AccuracyEvent & ek1=Interpol(trajecto.AccuracySeries().Event(Acc_evt1), trajecto.AccuracySeries().Event(Acc_evt1+1) ,time1) ; // Incertitude au temps t+1   evenement k+1

//        const double & dt=time1 - time;

//         Lg::Point3d Uncertainty_pos_ek=Propagation_uncertainty_from_NEB_to_LAM93_Position(m_Event.Event(k),ek);
//         Lg::Point3d Uncertainty_velocity_ek=Propagation_uncertainty_from_NEB_to_LAM93_Vitesse(m_Event.Event(k),ek);
//         Lg::Point3d Uncertainty_pos_ek1=Propagation_uncertainty_from_NEB_to_LAM93_Position(m_Event.Event(k+1),ek1); // Incertitude de vitesse de l'evenement k


//        v_normaliseX.push_back( facteur_normalisation(amplification_x*Uncertainty_pos_ek.x(),amplification_x*Uncertainty_pos_ek1.x(), Uncertainty_velocity_ek.x(),dt) );
//        v_normaliseY.push_back( facteur_normalisation(amplification_y*Uncertainty_pos_ek.y(),amplification_y*Uncertainty_pos_ek1.y(), Uncertainty_velocity_ek.y(),dt) );
////        v_normaliseZ.push_back( facteur_normalisation(amplification_z*Uncertainty_pos_ek.z(),amplification_z*Uncertainty_pos_ek1.z(), Uncertainty_velocity_ek.z(),dt) );

//    }

//       m_covarianceX=Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement);
//       m_covarianceY=Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement);
////       m_covarianceZ=Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement);

//    for(unsigned int i=0 ; i<m_NbreEvenement; ++i){

//        auto const& event_i = m_Event.Event(i);
//        const double & temps=event_i.m_time;
//        const int &  Accuracy_evt= trajecto.AccuracySeries().PrevIndex(temps);
//        auto const& ek=Interpol(trajecto.AccuracySeries().Event(Accuracy_evt), trajecto.AccuracySeries().Event(Accuracy_evt) ,temps) ; // Récupère l'incertitude au temps t

//        Lg::Point3d Uncertainty_pos_event=Propagation_uncertainty_from_NEB_to_LAM93_Position(m_Event.Event(i),ek);

//        const double & cxi = carre(amplification_x*Uncertainty_pos_event.x());
//        const double & cyi = carre(amplification_y*Uncertainty_pos_event.y());
////        const double czi = carre(amplification_z*Uncertainty_pos_event.z());

//        m_covarianceX(i,i) = cxi;
//        m_covarianceY(i,i) = cyi;
////        m_covarianceZ(i,i) = czi;

//        double normX =1.;
//        double normY = 1.;
////        double normZ= 1.;

//        for(unsigned int j=i+1 ; j<m_NbreEvenement; ++j){
//              // Construction matrice m_covarianceX
//            normX *= v_normaliseX[j-1];
//            m_covarianceX(i,j) = cxi *normX;
//            m_covarianceX(j,i) = m_covarianceX(i,j);
//              // Construction matrice m_covarianceY
//            normY *= v_normaliseY[j-1];
//            m_covarianceY(i,j) = cyi *normY ;
//            m_covarianceY(j,i) = m_covarianceY(i,j);

//             // Construction matrice m_covarianceZ
////            normZ*=v_normaliseZ[j-1];
////            m_covarianceZ(i,j) = czi *normZ  ;
////            m_covarianceZ(j,i)=m_covarianceZ(i,j);

//        }
//    }

//    v_normaliseX.clear();
//    v_normaliseY.clear();
////    v_normaliseZ.clear();

//    m_covarianceXY.resize(m_NbreEvenement*2,m_NbreEvenement*2);
//    m_covarianceXY.setZero();

////     Construction matrice m_covariance et m_covarianceXY
//     for(unsigned int i=0 ; i<m_NbreEvenement; ++i){
//           for(unsigned int j=i; j<m_NbreEvenement; ++j){

//               const double & Cx_ij=m_covarianceX(i,j);
//               const double & Cy_ij=m_covarianceY(i,j);

////               const double & Cz_ij=m_covarianceZ(i,j);
////               const unsigned int k=i*3;
////               const unsigned int p=j*3;

//               const unsigned int & I=2*i ;
//               const unsigned int & J=2*j ;

//               if(float(Cx_ij)!=0)
//                    m_covarianceXY.coeffRef(I,J) = Cx_ij;
//               if(float(Cy_ij)!=0)
//                    m_covarianceXY.coeffRef(I+1,J+1) = Cy_ij;

////               m_covariance(k,p)=Cx_ij;
////               m_covariance(k+1,p+1)=Cy_ij;
////               m_covariance(k+2,p+2)=Cz_ij;

//               if(i!=j){

//                   if(float(Cx_ij)!=0)
//                        m_covarianceXY.coeffRef(J,I) = Cx_ij;
//                   if(float(Cy_ij)!=0)
//                        m_covarianceXY.coeffRef(J+1,I+1) = Cy_ij;

////                   m_covariance(p,k)=Cx_ij;
////                   m_covariance(p+1,k+1)=Cy_ij;
////                   m_covariance(p+2,k+2)=Cz_ij;
//               }


//           }
//    }

//      std::cout << "Time covariance constructor " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

//}



//             // Construction matrice m_covarianceZ
////            normZ*=v_normaliseZ[j-1];
////            m_covarianceZ(i,j) = czi *normZ  ;
////            m_covarianceZ(j,i)=m_covarianceZ(i,j);

//        }
//    }

//    v_normaliseX.clear();
//    v_normaliseY.clear();
////    v_normaliseZ.clear();

//    m_covarianceXY.resize(m_NbreEvenement*2,m_NbreEvenement*2);
//    m_covarianceXY.setZero();

////     Construction matrice m_covariance et m_covarianceXY
//     for(unsigned int i=0 ; i<m_NbreEvenement; ++i){
//           for(unsigned int j=i; j<m_NbreEvenement; ++j){

//               const double & Cx_ij=m_covarianceX(i,j);
//               const double & Cy_ij=m_covarianceY(i,j);

////               const double & Cz_ij=m_covarianceZ(i,j);
////               const unsigned int k=i*3;
////               const unsigned int p=j*3;

//               const unsigned int & I=2*i ;
//               const unsigned int & J=2*j ;

//               if(float(Cx_ij)!=0)
//                    m_covarianceXY.coeffRef(I,J) = Cx_ij;
//               if(float(Cy_ij)!=0)
//                    m_covarianceXY.coeffRef(I+1,J+1) = Cy_ij;

////               m_covariance(k,p)=Cx_ij;
////               m_covariance(k+1,p+1)=Cy_ij;
////               m_covariance(k+2,p+2)=Cz_ij;

//               if(i!=j){

//                   if(float(Cx_ij)!=0)
//                        m_covarianceXY.coeffRef(J,I) = Cx_ij;
//                   if(float(Cy_ij)!=0)
//                        m_covarianceXY.coeffRef(J+1,I+1) = Cy_ij;

////                   m_covariance(p,k)=Cx_ij;
////                   m_covariance(p+1,k+1)=Cy_ij;
////                   m_covariance(p+2,k+2)=Cz_ij;
//               }


//           }
//    }

//      std::cout << "Time covariance constructor " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

//}


