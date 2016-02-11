#include "CovarianceMatrix.h"

template <typename T> inline constexpr T carre(T const& v) { return v * v; }

Covariance::Covariance(Trajecto &  trajecto) :
          m_AccEvent(trajecto.AccuracySeries()),
          m_NbreEvenement( m_AccEvent.Nevent() ),
          m_covarianceX(Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement)),
          m_covarianceY(Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement)),
          m_covarianceZ(Eigen::MatrixXd::Zero(m_NbreEvenement, m_NbreEvenement))
{
    std::vector<double> v_normaliseX;
    std::vector<double> v_normaliseY;
    std::vector<double> v_normaliseZ;

    assert(m_NbreEvenement > 0);
    for(int k=0 ; k<(m_NbreEvenement-1) ; ++k){

        auto const& ek1 = m_AccEvent.Event(k+1);
        auto const& ek = m_AccEvent.Event(k);

        const double dt=ek1.m_time - ek.m_time;

        v_normaliseX.push_back(facteur_normalisation(ek.m_northPositionRMSError,ek1.m_northPositionRMSError,ek.m_northVelocityRMSError,dt));
        v_normaliseY.push_back(facteur_normalisation(ek.m_eastPositionRMSError,ek1.m_eastPositionRMSError,ek.m_eastVelocityRMSError,dt));
        v_normaliseZ.push_back(facteur_normalisation(ek.m_downPositionRMSError,ek1.m_downPositionRMSError,ek.m_downVelocityRMSError,dt));
    }

    for(int i=0 ; i<m_NbreEvenement; ++i){

        auto const& event_i = m_AccEvent.Event(i);

        const double cxi = carre(event_i.m_northPositionRMSError);
        const double cyi = carre(event_i.m_eastPositionRMSError);
        const double czi = carre(event_i.m_downPositionRMSError);

        m_covarianceX(i,i) = cxi;
        m_covarianceY(i,i) = cyi;
        m_covarianceZ(i,i) = czi;

        double normX =1.;
        double normY = 1.;
        double normZ= 1.;

        for(int j=i+1 ; j<m_NbreEvenement; ++j){

            normX *= v_normaliseX[j-1];
            m_covarianceX(i,j) = cxi *normX;
            m_covarianceX(j,i) = m_covarianceX(i,j);

            normY *= v_normaliseY[j-1];
            m_covarianceY(i,j) = cyi *normY ;
            m_covarianceY(j,i) = m_covarianceY(i,j);

            normZ*=v_normaliseZ[j-1];
            m_covarianceZ(i,j) = czi *normZ  ;
            m_covarianceZ(j,i)=m_covarianceZ(i,j);

        }
    }
}






