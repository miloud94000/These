#ifndef __MODELBATI_H__
#define __MODELBATI_H__

//-----------------------------------------------------------------------------------------------------------//
/*!                                            Modele Bati3D!                                                */
//-----------------------------------------------------------------------------------------------------------//

#include <iostream>
#include <fstream>

#include <lib3ds/file.h>
#include <lib3ds/mesh.h>
#include <lib3ds/node.h>


#include "LgPoint3.hpp"
#include "LgPoint2.hpp"

#include "ClassKdTree.h"
#include "accelerators/kdtreeaccel.h"


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
#include <map>

using namespace std;

struct MyFace{
    uint16_t a, b, c;
};


struct strip{

    vector<Lg::Point3d> v_pt; // On stocke dans l'ordre les points A,B,C,D d'un patch
    unsigned int numero_strip;
};

struct MyMesh
{
    string nameBloc;
    string name;
    string nameBC;
    string nameBS;
    vector<Lg::Point3> pt;
    vector<Lg::Point2> uv;
    vector<MyFace> face;
    uint16_t first_id;
    vector<strip> v_strip;

};


class Modele_Bati3D
{


public :

   std::vector<Reference<Primitive>>  m_vref_triangle;       // Contient tous les triangles des facades et des toits (Format pour Ray-tracing)
   std::vector<KdTreeTriangleP*> m_vp_triangle;
   std::vector<Reference<Primitive>> m_vref_triangle_patch;  //Contient tous les triangles des strips et des toits  (Format pour Ray-tracing)
   std::vector<KdTreeTriangleP*> m_vp_triangle_patch;

   std::map<std::string,Lg::Point3d> m_map_facade_normale;





   Lg::Point3d m_pivot_BATI3D;                               // Point pivot du fichier 3DS
   double m_largeur_strip;

   Modele_Bati3D(Lib3dsFile *file,const Lg::Point3d & pivot_bati3D=Lg::Point3d(), const double & largeur_strip=1.5);
   MyMesh CreateMesh(Lib3dsMesh *mesh, string nameBC,string nameBS, string nameBloc);
   void MeshNode(Lib3dsFile *f, Lib3dsNode *node, vector<MyMesh> & v_mesh);


private :

    unsigned int m_numero_bande;
    vector<MyMesh> v_mesh;

};


#endif



#include "ModeleBati3D.h"

Modele_Bati3D::Modele_Bati3D(Lib3dsFile *file,const Lg::Point3d & pivot_bati3D, const double & largeur_strip) : m_pivot_BATI3D(pivot_bati3D), m_largeur_strip(largeur_strip) {

    m_numero_bande=1; // Chaque strip vas recevoir un numéro unique

    if (file==NULL )
    {

      std::cout<< RED <<__FUNCTION__<<" ....  "<< "Error loading: " << " => skipping" <<  "...."<<std::endl;

    }

    std::cout << "Loaded "<<std::endl;

    lib3ds_file_eval(file,0.0f); //time passed as second parameter (required to recompute coordinates)

    if (file->nodes == NULL)
    {

             std::cout<< RED <<__FUNCTION__<<" ....  "<< "Pas de noeuds" <<  "...."<<std::endl;
    }
    else
    {

            for(Lib3dsNode *node=file->nodes; node; node=node->next){
                MeshNode(file, node,v_mesh);
            }
    }


    Lg::Point3d normale_facade;

    for(unsigned i=0;i<v_mesh.size();i++){

        for(unsigned j=0;j<v_mesh[i].face.size();j++)
         {
              Lg::Point3d point_A=v_mesh[i].pt[v_mesh[i].face[j].a]+m_pivot_BATI3D;
              Lg::Point3d point_B=v_mesh[i].pt[v_mesh[i].face[j].b]+m_pivot_BATI3D;
              Lg::Point3d point_C=v_mesh[i].pt[v_mesh[i].face[j].c]+m_pivot_BATI3D;

              KdTreeTriangleP* triangle=new KdTreeTriangleP(point_A,point_B,point_C,v_mesh[i].nameBloc,v_mesh[i].nameBC,v_mesh[i].nameBS,v_mesh[i].name);
              normale_facade=triangle->GetTriangleNormal();
              m_vp_triangle.push_back(triangle);

         }

        m_map_facade_normale[v_mesh[i].name]=normale_facade;
    }

//    std::map<std::string,Lg::Point3d>::iterator it = m_map_facade_normale.begin();
//    for (it=m_map_facade_normale.begin(); it!=m_map_facade_normale.end(); ++it)
//       std::cout << it->first << " => " << it->second << '\n';



    CreateKdTreeTriangleVector(m_vp_triangle,m_vref_triangle);

    for(unsigned int i=0;i<v_mesh.size();i++){

        if(v_mesh[i].name[0]=='F'){
            for(unsigned int j=0;j<v_mesh[i].v_strip.size();j++)
             {

                      Lg::Point3d point_A_patch=v_mesh[i].v_strip[j].v_pt[0]+m_pivot_BATI3D;
                      Lg::Point3d point_B_patch=v_mesh[i].v_strip[j].v_pt[1]+m_pivot_BATI3D;
                      Lg::Point3d point_C_patch=v_mesh[i].v_strip[j].v_pt[2]+m_pivot_BATI3D;
                      Lg::Point3d point_D_patch=v_mesh[i].v_strip[j].v_pt[3]+m_pivot_BATI3D;

                     // Attention : les points doivent être entré dans le sens inverse des aiguilles d'une montre.

                      // A                B
                      // *---------------*
                      // |    Facade     |
                      // |               |
                      // *---------------*
                      // D                C

                      KdTreeTriangleP* triangle1=new KdTreeTriangleP(point_A_patch,point_C_patch,point_B_patch,v_mesh[i].nameBloc,v_mesh[i].nameBC,v_mesh[i].nameBS,v_mesh[i].name);
                      KdTreeTriangleP* triangle2=new KdTreeTriangleP(point_A_patch,point_D_patch,point_C_patch,v_mesh[i].nameBloc,v_mesh[i].nameBC,v_mesh[i].nameBS,v_mesh[i].name);

                      triangle1->setPatch(v_mesh[i].v_strip[j].numero_strip); // Chaque triangle à un numéro de strip
                      triangle2->setPatch(v_mesh[i].v_strip[j].numero_strip);

                      m_vp_triangle_patch.push_back(triangle1); //  : 1 strip = 2 triangles
                      m_vp_triangle_patch.push_back(triangle2);


              }

        }else{

            for(unsigned j=0;j<v_mesh[i].face.size();j++)
             {
                  Lg::Point3d point_A=v_mesh[i].pt[v_mesh[i].face[j].a]+m_pivot_BATI3D;
                  Lg::Point3d point_B=v_mesh[i].pt[v_mesh[i].face[j].b]+m_pivot_BATI3D;
                  Lg::Point3d point_C=v_mesh[i].pt[v_mesh[i].face[j].c]+m_pivot_BATI3D;


                    KdTreeTriangleP* triangle=new KdTreeTriangleP(point_A,point_B,point_C,v_mesh[i].nameBloc,v_mesh[i].nameBC,v_mesh[i].nameBS,v_mesh[i].name);


                   m_vp_triangle_patch.push_back(triangle);

                   point_A=point_A+pivot_bati3D;
                   point_B=point_B+pivot_bati3D;
                   point_C=point_C+pivot_bati3D;




             }




        }


    }

    v_mesh.clear();
    CreateKdTreeTriangleVector(m_vp_triangle_patch,m_vref_triangle_patch);
    std::cout<< GREEN <<__FUNCTION__<<"---Succès du chargement des données ----"<<DEFAULT_COLOR<<std::endl;

}


void Modele_Bati3D::MeshNode(Lib3dsFile *f, Lib3dsNode *node, vector<MyMesh> & v_mesh)
{

       Lib3dsNode *p;

       for (p=node->childs; p!=0; p=p->next){
           MeshNode(f, p, v_mesh);
       }

       Lib3dsMesh *mesh=lib3ds_file_mesh_by_name(f,node->name);

       if(!mesh) return;

       if(mesh->name[0]=='F' || mesh->name[0]=='T' ){
           string nameBloc=node->parent->parent->parent->parent->name;
           string nameBC=node->parent->parent->parent->name;
           string nameBS=node->parent->parent->name;

           if(nameBloc[0]=='B'){ /// Ajout de la condition car le fichier 3ds a parfois un défault de structure
                 v_mesh.push_back(CreateMesh(mesh,nameBC,nameBS,nameBloc));
           }else{
               cout << mesh->name <<endl;
           }

       }
}


MyMesh Modele_Bati3D::CreateMesh(Lib3dsMesh *mesh, string nameBC,string nameBS, string nameBloc){

    MyMesh my_mesh;

    std::string name(mesh->name);

    my_mesh.nameBloc = nameBloc;
    my_mesh.name = name;
    my_mesh.nameBC=nameBC;
    my_mesh.nameBS=nameBS;

    for (unsigned i=0; i<mesh->faces; ++i)
    {
        MyFace f;
        f.a = mesh->faceL[i].points[0];
        f.b = mesh->faceL[i].points[1];
        f.c = mesh->faceL[i].points[2];

        my_mesh.face.push_back(f);
    }


    for (unsigned i=0; i<mesh->points; ++i)
    {
        Lg::Point3d p;
        p.x() = mesh->pointL[i].pos[0];
        p.y() = mesh->pointL[i].pos[1];
        p.z() = mesh->pointL[i].pos[2];

        my_mesh.pt.push_back(p);
    }


    if(name[0]=='F'){

               Lg::Point3d A;
               Lg::Point3d B;
               Lg::Point3d C;
               Lg::Point3d D;

               // Dans le cas ou la facade à 4 ou 6 points
               if(my_mesh.pt.size()==4){  // Si la facade à 4 points

                   A=my_mesh.pt[0];
                   B=my_mesh.pt[1];
                   C=my_mesh.pt[2];
                   D=my_mesh.pt[3];

               }else if(my_mesh.pt.size()==6){ // Si la facade a 6 points



                   // A(pt0)    (pt1)       B (pt2)
                   // *-----     *    ----------*
                   // |    Facade               |
                   // |                         |
                   // *(pt5)-----*(pt4)   ---  -* (pt3)
                   // D                         C

                    A=my_mesh.pt[0];
                    B=my_mesh.pt[2];
                    C=my_mesh.pt[3];
                    D=my_mesh.pt[5];

               }else{
                       std::cout<< RED <<__FUNCTION__<<"WARNING : Les facades peuvent avoir plus de 6 points : revoir l'algo sur les patchs "<<DEFAULT_COLOR<<std::endl;

                       exit(1);
               }


                    // A                B
                    // *---------------*
                    // |    Facade     |
                    // |               |
                    // *---------------*
                    // D                C

                    const Lg::Point3d AB=B-A; // Vecteur de la longeur du haut de la facade
                    const Lg::Point3d DC=C-D; // Vecteur de la longeur du bas de la facade

                    const double Longueur_AB=AB.Norm();
                    const double Longueur_DC=DC.Norm();


                    const unsigned int nb_strip=round(Longueur_AB/m_largeur_strip);

                    const double fcd_strip_widthAB=Longueur_AB/nb_strip;
                    const double fcd_strip_widthDC=Longueur_DC/nb_strip;


                    const Lg::Point3 AB_unitaire=AB.Normalized();
                    const Lg::Point3 DC_unitaire=DC.Normalized();

                    if(nb_strip==0){

                       strip bande;
                       bande.v_pt.push_back(A);
                       bande.v_pt.push_back(B);
                       bande.v_pt.push_back(C);
                       bande.v_pt.push_back(D);

                       bande.numero_strip=m_numero_bande;
                       m_numero_bande++;

                       my_mesh.v_strip.push_back(bande);

                    }else{


                        Lg::Point3 a,b,c,d; // Les 4 points de la bande

                        for(unsigned int i=1;i<=nb_strip;i++){
                                strip bande;

                                if(i==1){

                                    a = A;
                                    b = A + AB_unitaire * fcd_strip_widthAB;
                                    d= D;
                                    c= D + DC_unitaire * fcd_strip_widthDC;

                                }else if(i==nb_strip){

                                    a=b;
                                    b=B;
                                    d=c;
                                    c=C;


                               }else{


                                 a=b;
                                 b= a + AB_unitaire * fcd_strip_widthAB;
                                 d=c;
                                 c=d +DC_unitaire * fcd_strip_widthDC;

                               }

                                bande.v_pt.push_back(a);
                                bande.v_pt.push_back(b);
                                bande.v_pt.push_back(c);
                                bande.v_pt.push_back(d);

                                bande.numero_strip=m_numero_bande;
                                m_numero_bande++;
                                my_mesh.v_strip.push_back(bande);

                     }
                 }



    }




    return my_mesh;

}


