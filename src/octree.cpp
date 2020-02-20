#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Vertex_handle Vertex_handle;


double max_sommets = 10;

class Noeud{
    private: 
        std::vector<Point_3> boite;
        std::map<Vertex_handle, Point_3> sommets;
        std::vector<Point_3> coord_noeud;

    public:
        Noeud(std::vector<Point_3> coord_noeud_mere) {
            for(int i = 0; i < coord_noeud_mere.size(); i++) {
                coord_noeud.push_back(coord_noeud_mere[i]);
            }
        }

        void setPoint(double x, double y, double z){
            Point_3 point(x,y,z);
            boite.push_back(point);
        }

        void setBoite(double min_x,double max_x,double min_y,double max_y,double min_z,double max_z){
            setPoint(min_x, min_y, min_z);
            setPoint(min_x, max_y, min_z);
            setPoint(min_x, min_y, max_z);
            setPoint(min_x, max_y, max_z);
            setPoint(max_x, min_y, min_z);
            setPoint(max_x, max_y, min_z);
            setPoint(max_x, max_y, max_z);
            setPoint(max_x, min_y, max_z);
        }

        void setCoord_noeud(Point_3 coords) {
            coord_noeud.push_back(coords);
        }

        void addSommet(Vertex_handle sommet, Point_3 coords) {
            sommets.insert(std::pair<Vertex_handle, Point_3>(sommet, coords));
        }


        std::vector<Point_3> getBoite() const{
            return boite;
        }
        
        std::map<Vertex_handle, Point_3> getSommets() const{
            return sommets;
        }

        std::vector<Point_3> getCoord() const{
            return coord_noeud;
        }


};

std::map<Vertex_handle, Point_3> calcul_sommets(Polyhedron & mesh) {
    std::map<Vertex_handle, Point_3> sommets;
    for(auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
        sommets[i] = i->point();
    }
    //Affiche les différents périmètres sur la sortie standard
    //for(auto it : sommets) std::cout << "numero sommet: " << it.first << " " << it.second[1] << "\n";
    return sommets;
}



void octree(std::map<Vertex_handle, Point_3> map_sommets, std::vector<Point_3> coord_noeud_mere) {

    int x , y, z;
    //std::map<Vertex_handle, Point_3> classification;
    if(map_sommets.size() > max_sommets && deep <6) {
        Noeud fils_1 = Noeud(coord_noeud_mere);
        Noeud fils_2 = Noeud(coord_noeud_mere);
        Noeud fils_3 = Noeud(coord_noeud_mere);
        Noeud fils_4 = Noeud(coord_noeud_mere);
        Noeud fils_5 = Noeud(coord_noeud_mere);
        Noeud fils_6 = Noeud(coord_noeud_mere);
        Noeud fils_7 = Noeud(coord_noeud_mere);
        Noeud fils_8 = Noeud(coord_noeud_mere);

        double min_x = map_sommets[0][0];
        double min_y = map_sommets[0][1];
        double min_z = map_sommets[0][2];
        double max_x = map_sommets[0][0];
        double max_y = map_sommets[0][1];
        double max_z = map_sommets[0][2];
        
        for(auto it : map_sommets) {
            if (it.second[0] > max_x) max_x = it.second[0];
            if (it.second[1] > max_y) max_y = it.second[1];
            if (it.second[2] > max_z) max_z = it.second[2];
            if (it.second[0] < min_x) min_x = it.second[0];
            if (it.second[1] < min_y) min_y = it.second[1];
            if (it.second[2] < min_z) min_z = it.second[2];
        }
    

        fils_1.setBoite(min_x, (max_x-min_x)/2, min_y, (max_y-min_y)/2, min_z, (max_z-min_z)/2);
        fils_2.setBoite(min_x, (max_x-min_x)/2, (max_y-min_y)/2, max_y, min_z, (max_z-min_z)/2);
        fils_3.setBoite(min_x, (max_x-min_x)/2, min_y, (max_y-min_y)/2, (max_z-min_z)/2, max_z);
        fils_4.setBoite(min_x, (max_x-min_x)/2, (max_y-min_y)/2, max_y, (max_z-min_z)/2, max_z);
        fils_5.setBoite((max_x-min_x)/2, max_x, min_y, (max_y-min_y)/2, min_z, (max_z-min_z)/2);
        fils_6.setBoite((max_x-min_x)/2, max_x, (max_y-min_y)/2, max_y, min_z, (max_z-min_z)/2);
        fils_7.setBoite((max_x-min_x)/2, max_x, min_y, (max_y-min_y)/2, (max_z-min_z)/2, max_z);
        fils_8.setBoite((max_x-min_x)/2, max_x, (max_y-min_y)/2, max_y, (max_z-min_z)/2, max_z);

        Point_3 point_central((max_x + min_x)/2.0,(max_y + min_y)/2.0,(max_z + min_z)/2.0);
        
        for(auto it : map_sommets) {
            x = 0;
            y = 0;
            z = 0;
            if(it.second[0] > point_central[0]) {
                x=  1;
            }
        
            if(it.second[1] > point_central[1]) {
                y = 1;
            }
         
            if(it.second[2] > point_central[2]) {
                z = 1;
            }
        }
    }
    else{
        std::cout<<"octree fini."<<std::endl;
    }
}

std::vector< Noeud > generateOctree(std::map<Vertex_handle, Point_3> map_sommets){
    std::vector< Noeud > address;
    octree(map_sommets,address,0);

    return address;
}


int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Il manque un parametre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
        return 1;
    }
   
    Polyhedron mesh;
    std::ifstream input(argv[1]);
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
        return 1;
    }
 
    unsigned int nbVerts = 0;
    for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
        ++nbVerts;
    }
    std::cout << "Nombre de sommets: " << nbVerts << std::endl;
   
    unsigned int nbEdges = 0;
    for (Halfedge_iterator i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i) {
        ++nbEdges;
    }
    nbEdges /= 2;
    std::cout << "Nombre d'aretes: " << nbEdges << std::endl;
 
    unsigned int nbFaces = 0;
    for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
        ++nbFaces;
    }
    std::cout << "Nombre de faces: " << nbFaces << std::endl;
   
    unsigned int euler = nbVerts - nbEdges + nbFaces;
    unsigned int genus = (2 - euler) / 2;
    std::cout << "En supposant que le maillage contienne une unique surface sans bord, alors son genre est de " << genus << std::endl;

    std::map<Vertex_handle, Point_3> map_sommets = calcul_sommets(mesh);
    //octree(map_sommets);

    return 0;
}

