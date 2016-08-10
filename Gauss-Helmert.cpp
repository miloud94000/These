#include "GaussHelmertModel.h"




Gauss_Helmert_Model::Gauss_Helmert_Model(const Ray_tracing & ray_tracing, XMls & mls, VarianceModelBati3D variance) :
nombre_de_bloc(ray_tracing.mv_bloc.size()),nombre_de_facade(ray_tracing.mv_facade.size()),
nombre_de_point(ray_tracing.m_vp_triangle.size()), nombre_de_patch(ray_tracing.mv_numero_strip.size()), nombre_event(mls.m_trajecto.AccuracySeries().Nevent() /*on travaille à 1Hz*/),
sigmaFacade(variance.sigmaFacade),sigmaBloc(variance.sigmaBloc),sigmaStrip(variance.sigmaStrip)
{


    cout << "Nombre de point :" <<nombre_de_point<<endl;
    setB(ray_tracing,mls); // Initialise la matrice m_B
    setCll(ray_tracing,mls); // Initialise la matrice m_Cll
    set_gl0(ray_tracing); // Initialise la fonction m_gl0
    set_l();
    set_l0(); //Initialise la fonction m_l et m_l0


   m_Cg=-1*m_gl0;

   Eigen::SparseMatrix<double> M = m_B_transpose*m_Cll*m_B;
{
//   Eigen::MatrixXd M_dense=toDense(M);

//    Eigen::MatrixXd M_dense_inverse=M_dense.inverse();
//    cout <<"Determiant dense" <<M_dense.determinant()<<endl;

//    Eigen::MatrixXd identite( M_dense.rows() , M_dense.rows() );

//    identite.setIdentity();

//   Eigen::MatrixXd AA=(M_dense*M_dense_inverse-identite);


//   cout << "Max coeff AA"<<AA.maxCoeff() <<endl;
//   cout << "min coeff AA"<<AA.minCoeff() <<endl;
}


   Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
   solver.compute(M);


   if(solver.info()!=Eigen::Success) {
        cout << "Gauss-Helmert Non resolvable" << endl;
   }

   cout << "Determinant : "<<solver.determinant() << endl;

   Eigen::SparseMatrix<double> Id(M.rows(),M.rows());
   Id.setIdentity();

    M_inverse = solver.solve(Id);



    m_v=m_Cll*m_B*M_inverse*m_Cg;




    m_l=m_l0+m_v;
    m_Cg=-m_gl0-(m_B_transpose*m_v);
    m_v=m_Cll*m_B*(M_inverse)*m_Cg;

    m_l=m_l+m_v;


     cout << m_v.maxCoeff()<<endl;



        //TODO : Faire une classe recalage
      cout <<m_l <<endl;
//    cout <<"---------------------------"<<endl;
//    cout <<"---------------------------"<<endl;

    for(unsigned int i=0; i<nombre_event;i++){

            Lg::Point2d dxdy(m_l[2*i],m_l[2*i+1]);
            v_Delta_trajectoire_xy[i]=dxdy;
           cout << "Trajectoire  :  " <<  i << dxdy <<endl;

    }



    for(unsigned int i=0; i<nombre_de_facade;i++){

            cout <<ray_tracing.mv_facade[i]<<" : "<< m_l[nombre_event*2+i] <<endl;
            v_offset_facade[ray_tracing.mv_facade[i]]=m_l[nombre_event*2+i];

    }

    cout <<"---------------------------"<<endl;
    cout <<"---------------------------"<<endl;
    for(unsigned int i=0; i<nombre_de_bloc;i++){

            Lg::Point2d dxdy(m_l[nombre_event*2+nombre_de_facade+2*i],m_l[nombre_event*2+nombre_de_facade+2*i+1]);
            v_offset_bloc_xy[ray_tracing.mv_bloc[i]]=dxdy;
            cout <<ray_tracing.mv_bloc[i]<<" : "<< m_l[nombre_event*2+nombre_de_facade+2*i] <<", "<< m_l[nombre_event*2+nombre_de_facade+2*i+1]<<endl;

    }


    effacer_matrices(); // Libère la mémoire quand nous avons trouvé la solution


}


void Gauss_Helmert_Model::setB(const Ray_tracing & ray_tracing,XMls & mls){


      Eigen::SparseMatrix<double> B(nombre_de_patch,nombre_event*2); // Dérivée de g par rapport à la trajectoire à 100Hz
      B.setZero();

      Eigen::SparseMatrix<double> C(nombre_de_patch,nombre_de_facade); // Dérivée de g par rapport au facade
      C.setZero();
      Eigen::SparseMatrix<double> D(nombre_de_patch,nombre_de_bloc*2); // Dérivée de g par rapport au bloc  btx bty
      D.setZero();
      Eigen::SparseMatrix<double> E(nombre_de_patch,nombre_de_patch);  // Dérivée de g par rapport au patch
      E.setZero();


      for(unsigned int k =0; k<nombre_de_point;k++){

           KdTreeTriangleP * KdTreeTrianglePoint = ray_tracing.m_vp_triangle[k]; // L'appariemment point-triangle issue du ray-tracing

           const unsigned int & j = index_strip(ray_tracing.mv_numero_strip, KdTreeTrianglePoint->numero_patch);         // Indice du strip
           const unsigned int & i = KdTreeTrianglePoint->getPrevIndex();                                                // Indice de l'évenement T- de la trajectoire

           const double & alpha = KdTreeTrianglePoint->alpha;                                                 //  Un point P= Pi + (1-alpha) delta- + alpha delta+
           const Lg::Point3d & nf = KdTreeTrianglePoint->GetTriangleNormal();

           Lg::Point3d derivation_by_trajectory_i = (alpha-1)*nf;  // Dérivée de g par rapport à T-
           Lg::Point3d derivation_by_trajectory_i1 = (-alpha)*nf;   // Dérivée de g par rapport à T+

           B.coeffRef(j,2*i)   += derivation_by_trajectory_i.x();
           B.coeffRef(j,2*i+1) += derivation_by_trajectory_i.y();
           B.coeffRef(j,2*i+2) += derivation_by_trajectory_i1.x();
           B.coeffRef(j,2*i+3) += derivation_by_trajectory_i1.y();

           const unsigned int & b = index_in_vector_string(ray_tracing.mv_facade, KdTreeTrianglePoint->m_nom_facade);    // Indice de la facade ray_tracing.mv_facade-> Les facades sont dans un certain ordre
           C.coeffRef(j,b) += 1.;

           const unsigned int & d = index_in_vector_string(ray_tracing.mv_bloc, KdTreeTrianglePoint->m_nom_du_bloc);
           D.coeffRef(j,2*d) += nf.x();
           D.coeffRef(j,2*d+1) += nf.y();


           E.coeffRef(j,j) += 1.;
        }

        Eigen::SparseMatrix<double> BCDE=hcat(B,C,D,E);

        m_B_transpose=BCDE;

        m_B= BCDE.transpose();
}

void Gauss_Helmert_Model::setCll(const Ray_tracing & ray_tracing , XMls & mls){

       Covariance covar(mls.m_trajecto,3.,3.,3.);

       Eigen::MatrixXd Cll_trajectoire_dense=covar.getCovarianceXY();

       Eigen::SparseMatrix<double> Cll_trajectoire=toSparse(Cll_trajectoire_dense);


       Eigen::SparseMatrix<double> Cll_facade(nombre_de_facade,nombre_de_facade);
       Cll_facade.setZero();

       for(unsigned int i =0; i<nombre_de_facade;i++){

          Cll_facade.coeffRef(i,i)=sigmaFacade*sigmaFacade;

       }

       Eigen::SparseMatrix<double> Cll_bloc(nombre_de_bloc*2,nombre_de_bloc*2);
       Cll_bloc.setZero();

       for(unsigned int i =0; i<nombre_de_bloc;i++){

           const unsigned int & k=2*i;
           Cll_bloc.coeffRef(k,k)=sigmaBloc*sigmaBloc;
           Cll_bloc.coeffRef(k+1,k+1)=sigmaBloc*sigmaBloc;

       }

       Eigen::SparseMatrix<double> Cll_strip(nombre_de_patch,nombre_de_patch);
       Cll_strip.setZero();

       for(unsigned int i =0; i<nombre_de_patch;i++){

             Cll_strip.coeffRef(i,i)=sigmaStrip*sigmaStrip;

        }

        m_Cll=blkdiag(Cll_trajectoire,Cll_facade, Cll_bloc,Cll_strip);; // Prend les 4 matrices et construit une nouvelle matrice avec à la diagonale ces 4 matrices




}

void Gauss_Helmert_Model::set_gl0(const Ray_tracing & ray_tracing){

     Eigen::VectorXd g_l0(nombre_de_patch);
     g_l0.setZero();

     for(unsigned int i =0; i<nombre_de_point;i++){

        KdTreeTriangleP * trianglePoint = ray_tracing.m_vp_triangle[i];
        Lg::Point3d Point_facade = trianglePoint->A();                  // Correspond à un point du strip
        Lg::Point3d Point_laser = trianglePoint->Pi;                    // Correspond au point apparié au strip
        Lg::Point3d normale_patch = trianglePoint->GetTriangleNormal(); //  Normale au strip

        const unsigned int & numero_patch = trianglePoint->numero_patch;
        const unsigned int & j = index_strip(ray_tracing.mv_numero_strip ,numero_patch );  // Position dans le vecteur

        g_l0.coeffRef(j) += (Point_facade-Point_laser)*normale_patch;

     }

     m_gl0=g_l0;


}

void Gauss_Helmert_Model::set_l()
{
     Eigen::VectorXd l(nombre_event*2+nombre_de_facade+nombre_de_bloc*2+nombre_de_patch);
     l.setZero();

     m_l=l;
}


void Gauss_Helmert_Model::set_l0()
{
    Eigen::VectorXd l0(nombre_event*2+nombre_de_facade+nombre_de_bloc*2+nombre_de_patch);
    l0.setZero();

    m_l0=l0;
}

void Gauss_Helmert_Model::effacer_matrices(){

    m_gl0.resize(0,0);
    m_B.resize(0,0);
    m_B_transpose.resize(0,0);
    m_Cll.resize(0,0);
    M_inverse.resize(0,0);
    m_Cg.resize(0,0);
    m_v.resize(0,0);


}






