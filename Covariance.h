#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <functional>
#include <numeric>
#include "ReaderFilePosPac/trajecto_reader.h"


/*! La class CovarianceMatrix a pour but de calculer les matrices de covariances X,Y et Z à partir d'une trajectoire */

class CovarianceMatrix
{


public:


    CovarianceMatrix(const Trajecto & trajecto) :
      m_trajecto{trajecto}
    {

        Covariance();
    }

    Eigen::MatrixXd getCovarianceX(){return m_covarianceX;}
    Eigen::MatrixXd getCovarianceY(){return m_covarianceY;}
    Eigen::MatrixXd getCovarianceZ(){return m_covarianceZ;}


private :

    /*! Objet Trajectoire contenant toutes les informations (Temps t, incertitudes positions et vitesses) pour calculer les matrices de covariances */
    Trajecto m_trajecto;

    /*! Matrice de covariance calculé à partir du temps t et des incertitudes positions, vitesses dans la direction X */
    Eigen::MatrixXd m_covarianceX;

    /*! Matrice de covariance calculé à partir du temps t  et des incertitudes positions, vitesses dans la direction Y */
    Eigen::MatrixXd m_covarianceY;

    /*! Matrice de covariance calculé à partir du temps t et  des incertitudes positions, vitesses dans la direction Z*/
    Eigen::MatrixXd m_covarianceZ;


     void Covariance(); /*! Permet de calculer les attributs m_covarianceX, m_covarianceY et m_covarianceZ*/



    inline double facteur_normalisation(const double & sigma_p_i,const double & sigma_p_i_1,const double  & sigma_v_i,const double & dt){

         return sqrt((sigma_p_i_1*sigma_p_i_1)/ ((dt*dt*sigma_v_i*sigma_v_i)+(sigma_p_i*sigma_p_i)) );
    }







};
