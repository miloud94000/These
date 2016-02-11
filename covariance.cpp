#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <functional>
#include <numeric>
#include "ReaderFilePosPac/trajecto_reader.h"


/*! La class CovarianceMatrix a pour but de calculer les matrices de covariances X,Y et Z à partir d'une trajectoire */

class Covariance
{


public:


    Covariance(Trajecto & trajecto);


     Eigen::MatrixXd const& getCovarianceX(){return m_covarianceX;}
     Eigen::MatrixXd const& getCovarianceY(){return m_covarianceY;}
     Eigen::MatrixXd const& getCovarianceZ(){return m_covarianceZ;}


private :

    /*! Objet Trajectoire contenant toutes les informations (Temps t, incertitudes positions et vitesses) pour calculer les matrices de covariances */



    TEventSeries<AccuracyEvent> & m_AccEvent;
    int m_NbreEvenement;
    Eigen::MatrixXd  m_covarianceX;
    Eigen::MatrixXd m_covarianceY;
    Eigen::MatrixXd m_covarianceZ;



    double facteur_normalisation(const double & sigma_p_i,const double & sigma_p_i_1,const double  & sigma_v_i,const double & dt){

         return sqrt((sigma_p_i_1*sigma_p_i_1)/ ((dt*dt*sigma_v_i*sigma_v_i)+(sigma_p_i*sigma_p_i)) );
    }







};
