#ifndef GAUSSHELMERT_H
#define GAUSSHELMERT_H

#include <Eigen/LU>
#include <iostream>
#include <fstream>
#include "LgPoint3.hpp"
#include "LgPoint2.hpp"
#include "ClassKdTree.h"
#include "XMlsUnc.h"
#include "CovarianceMatrix.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include "RayTracing.h"

#define COUT_SUB 1000000
#define BLACK    "\033[1;30m"
#define RED      "\033[1;31m"
#define GREEN    "\033[1;32m"
#define YELLOW   "\033[1;33m"
#define BLUE     "\033[1;34m"
#define PURPLE   "\033[1;35m"
#define CYAN     "\033[1;36m"
#define GREY     "\033[1;37m"
#define DEFAULT_COLOR "\033[0;m"



using namespace std;

struct VarianceModelBati3D{

    double sigmaFacade; // Incertitude identique pour toutes les facades (Necessaire dans la matrice de covariance des observations) Défault 10 cm
    double sigmaBloc;   // Incertitude identique pour tous les blocs (Necessaire dans la matrice de covariance des observations)      Défault 1m
    double sigmaStrip;  // Incertitude identique pour tous les Strip (Necessaire dans la matrice de covariance des observations)     Défault 20 cm


    VarianceModelBati3D(): sigmaFacade(0.05), sigmaBloc(0.1), sigmaStrip(0.05){}

};


class Gauss_Helmert_Model{

public :
    /*!

Model of constraint between the observations only : g(l)=0

     g(l)=∑(Q -P +T ).nf +of + s =0

     where :
     T_{b} block translation T_{b}=[Tx Ty Tz]^{T} and b is the block index, Tx and Ty are the observations
     Q^{b}_{j} Point of the facade
     P_{i} is the laser point The are three observations : The position of the vehicule (in the direction x, y and Z)
     t_j is an scalar. And it's an observation

     Linear Model : B^T v = Cg
     where:
           B  = ((∂g/ ∂l)|l=l(0))^T is the matrix of partial derivatives with respect to observations
           v  = Cll B (B^T Cll B)^(-1) Cg
           Cg = -g(l(0))- B^T (l-l(0)) is the correction vector
           Cll is the covariance matrice of observations
           l is the observations vector and l(0) is initial vector

    !*/

    // Solution au probleme de Gauss Helmert
    std::map<std::string,double> v_offset_facade; // Solution sur les observations des offsets de facades
    std::map<std::string,Lg::Point2d> v_offset_bloc_xy;  // Solution sur les observations des offsets de blocs
    std::map<unsigned int ,Lg::Point2d> v_Delta_trajectoire_xy;  // //Solution sur les observations des évenement de  la trajectoire



    Eigen::SparseMatrix<double>  m_Cll;// Matrice creuse

    Eigen::SparseMatrix<double> m_B_transpose,m_B;
    Eigen::VectorXd m_l,m_l0,m_Cg,m_gl0,m_v;
    Eigen::SparseMatrix<double> M_inverse;
    unsigned int nombre_de_bloc;
    unsigned int nombre_de_facade;
    unsigned int nombre_de_point;
    unsigned int nombre_de_patch;
    unsigned int nombre_event;


    double sigmaFacade; // Incertitude identique pour toutes les facades (Necessaire dans la matrice de covariance des observations) Défault 10 cm
    double sigmaBloc;   // Incertitude identique pour tous les blocs (Necessaire dans la matrice de covariance des observations)      Défault 1m
    double sigmaStrip;  // Incertitude identique pour tous les Strip (Necessaire dans la matrice de covariance des observations)     Défault 20 cm

    Gauss_Helmert_Model(const Ray_tracing & ray_tracing, XMls & mls, VarianceModelBati3D variance=VarianceModelBati3D());

private :

    void setB(const Ray_tracing & ray_tracing,XMls & mls);
    void setCll(const Ray_tracing & ray_tracing , XMls & mls);
    void set_gl0(const Ray_tracing & ray_tracing);
    void set_l();
    void set_l0();
    void effacer_matrices();

inline Eigen::MatrixXd toDense(const Eigen::SparseMatrix<double> &input) {
   Eigen::MatrixXd output(input.rows(), input.cols());
   for (int k = 0; k < input.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(input,k); it; ++it) {
          double value = input.coeff(it.row(), it.col());
          output(it.row(), it.col()) = value;
      }
    }
    return output;

}



      // Ou se trouve value (unsigned int) dans le vecteur ? Les strip sont rangés dans un certain ordre celle issue du ray_tracing
inline unsigned int index_strip(const vector<unsigned int> & mv_numero_strip,const unsigned int & value){

    auto o=find(mv_numero_strip.begin(),mv_numero_strip.end(),value);
    unsigned int indice_strip;
    if (o != mv_numero_strip.end()){

          indice_strip= distance(mv_numero_strip.begin(), o);

     }else{
                 std::cout<< RED <<__FUNCTION__<<" ....  "<< "Le patch "<< value << " n'est pas dans le vecteur (mv_numero_strip) "  <<  "...."<<std::endl;
                 exit(1);
     }

    return indice_strip ;
}

     // Ou se trouve value (string) dans le vecteur ?
inline unsigned int index_in_vector_string(const vector<std::string> & mv,const std::string & value){

     auto o1=find(mv.begin(),mv.end(),value);
     unsigned int indice;

     if (o1 != mv.end())
     {
         indice= distance(mv.begin(), o1); // On récupère l'index de la facade

      }else{

          std::cout<< RED <<__FUNCTION__<<" ....  "<< "La valeur "<< value << " n'est pas dans le vecteur (mv) "  <<  "...."<<std::endl;
          exit(1);
      }

     return indice;
}

     //Permet de concaténer 2 matrices creuses à l'horizontale
inline void hcat(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, Eigen::SparseMatrix<double>& result) {
      assert(A.rows() == B.rows());
      result = Eigen::SparseMatrix<double>(A.rows(), A.cols() + B.cols());
      result.middleCols(0,A.cols()) = A;
      result.middleCols(A.cols(), B.cols()) = B;
}

inline Eigen::SparseMatrix<double> hcat(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B) {

      Eigen::SparseMatrix<double> r;
      hcat(A, B, r);
      return r;
}

     // Concaténation horizontales des matrices creuses A B C D
inline Eigen::SparseMatrix<double> hcat(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B,const Eigen::SparseMatrix<double>& C,const Eigen::SparseMatrix<double>& D) {

      Eigen::SparseMatrix<double> r1,r2,r3;

      hcat(A, B, r1);    // r1=[A|B]
      hcat(C, D, r2);   // r2=[C|D]
      hcat(r1,r2,r3);   // r3=[r1|r2]=[A|B|C|D]

      return r3;
}

inline Eigen::SparseMatrix<double> toSparse(const Eigen::MatrixXd &input) {

   Eigen::SparseMatrix<double> output(input.rows(), input.cols());
   for(int i = 0; i < input.rows(); i++) {
       for(int j = 0; j < input.cols(); j++) {
            if (abs(input(i,j)) > 1e-30) {
              output.coeffRef(i,j) = input(i,j);
            }
       }
    }
    return output;
}

       //Permet de concaténer 2 matrices creuses à la verticale
inline void vcat(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, Eigen::SparseMatrix<double,Eigen::RowMajor>& result) {
    assert(A.cols() == B.cols());
    result = Eigen::SparseMatrix<double,Eigen::RowMajor>(A.rows() + B.rows(), A.cols());

    result.middleRows(0,A.rows()) = A;
    result.middleRows(A.rows(), B.rows()) = B;
}

inline Eigen::SparseMatrix<double> vcat(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B) {
    Eigen::SparseMatrix<double,Eigen::RowMajor> r;
    vcat(A, B, r);

    return r;
}

       //Permet de construire une matrice creuse  en mettant 2 matrices en block à la diagonale
inline void blkdiag(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, Eigen::SparseMatrix<double,Eigen::RowMajor>& result) {
       // concatenate A horizontally with a block of zeros of size (A.rows, B.cols)
   Eigen::SparseMatrix<double>block21(A.rows(), B.cols());
   Eigen::SparseMatrix<double>block1;
   hcat(A, block21, block1);
       // concatenate B horizontally with a block of zeros of size (B.rows, A.cols)
   Eigen::SparseMatrix<double>block12(B.rows(), A.cols());
   Eigen::SparseMatrix<double>block2;
   hcat(block12, B, block2);
       // vertically concatenate the two blocks
    vcat(block1, block2, result);
}

inline Eigen::SparseMatrix<double> blkdiag(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B) {

  Eigen::SparseMatrix<double,Eigen::RowMajor> r;
  blkdiag(A, B, r);
  return r;

}

inline Eigen::SparseMatrix<double> blkdiag(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B,const Eigen::SparseMatrix<double>& C,const Eigen::SparseMatrix<double>& D) {

  Eigen::SparseMatrix<double,Eigen::RowMajor> r1,r2,r3;
  blkdiag(A, B, r1);
  blkdiag(C, D, r2);
  blkdiag(r1,r2,r3);

  return r3;
}







};
#endif
