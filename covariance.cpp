#include "CovarianceMatrix.h"


void CovarianceMatrix::CovarianceMatrix(){

    unsigned NbreEvenement =m_trajecto.AccuracySeries().Nevent(); // Indique le nombre d'évenement de la trajectoire

    TEventSeries<AccuracyEvent> AccEvent=m_trajecto.AccuracySeries();

    std::vector<double> v_normaliseX; /* Facteur de normalisation suivant X */
    std::vector<double> v_normaliseY; /* Facteur de normalisation suivant Y */
    std::vector<double> v_normaliseZ; /* Facteur de normalisation suivant Z */

    double dt;

    for(unsigned i=0;i<(NbreEvenement-1);i++){

        dt=AccEvent.Event(i+1).m_time-AccEvent.Event(i).m_time; // Différence de temps entre 2 évenements

        /* Facteur de normalisation Ni suivant X */
        v_normaliseX.push_back(facteur_normalisation(AccEvent.Event(i).m_northPositionRMSError,AccEvent.Event(i+1).m_northPositionRMSError,AccEvent.Event(i).m_northVelocityRMSError,dt));

        /* Facteur de normalisation Ni suivant X */
        v_normaliseY.push_back(facteur_normalisation(AccEvent.Event(i).m_eastPositionRMSError,AccEvent.Event(i+1).m_eastPositionRMSError,AccEvent.Event(i).m_eastVelocityRMSError,dt));

        /* Facteur de normalisation Ni suivant Z */
        v_normaliseZ.push_back(facteur_normalisation(AccEvent.Event(i).m_downPositionRMSError,AccEvent.Event(i+1).m_downPositionRMSError,AccEvent.Event(i).m_downVelocityRMSError,dt));

    }

    m_covarianceX=Eigen::MatrixXd::Zero(NbreEvenement, NbreEvenement); /* Matrice carré */
    m_covarianceY=Eigen::MatrixXd::Zero(NbreEvenement, NbreEvenement); /* Matrice carré */
    m_covarianceZ=Eigen::MatrixXd::Zero(NbreEvenement, NbreEvenement); /* Matrice carré */


   double SommeNix,SommeNiy,SommeNiz ;

   for(unsigned i=0;i<NbreEvenement;i++){

       for(unsigned j=i;j<NbreEvenement;j++){

           /*si i==j alors cov(Di,Dj)=(Sigma_position_i)²*/
           m_covarianceX(i,j)=AccEvent.Event(i).m_northPositionRMSError*AccEvent.Event(i).m_northPositionRMSError; /*(Sigma_px_i)² */
           m_covarianceY(i,j)=AccEvent.Event(i).m_eastPositionRMSError*AccEvent.Event(i).m_eastPositionRMSError;   /*(Sigma_py_i)² */
           m_covarianceZ(i,j)=AccEvent.Event(i).m_downPositionRMSError*AccEvent.Event(i).m_downPositionRMSError;   /*(Sigma_pz_i)² */

           if(i!=j){

               /*si i!=j alors cov(Di,Dj)=(Sigma_position_i)² ∑ facteur de normalisation de i à j-1  */

               SommeNix=accumulate(v_normaliseX.begin()+i, v_normaliseX.begin()+j, 1., std::multiplies<double>());
               m_covarianceX(i,j)=m_covarianceX(i,j)*SommeNix;
               m_covarianceX(j,i)=m_covarianceX(i,j); /* Matrice symétrique */

               SommeNiy=accumulate(v_normaliseY.begin()+i, v_normaliseY.begin()+j, 1., std::multiplies<double>());
               m_covarianceY(i,j)=m_covarianceY(i,j)*SommeNiy;
               m_covarianceY(j,i)=m_covarianceY(i,j); /* Matrice symétrique */

               SommeNiz=accumulate(v_normaliseZ.begin()+i, v_normaliseZ.begin()+j, 1., std::multiplies<double>());
               m_covarianceZ(i,j)=m_covarianceZ(i,j)*SommeNiz;
               m_covarianceZ(j,i)=m_covarianceZ(i,j); /* Matrice symétrique */

           }


       }



   }



}










