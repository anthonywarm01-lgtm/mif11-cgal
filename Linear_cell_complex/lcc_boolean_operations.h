#ifndef LCC_BOOLEAN_OPERATIONS_H
#define LCC_BOOLEAN_OPERATIONS_H

#include<array>
#include<queue>
#include<map>
#include<unordered_map>
#include<vector>
#include<CGAL/Surface_mesh.h>
#include<CGAL/Kernel_traits.h>
#include<CGAL/Polygon_mesh_processing/autorefinement.h>
#include<CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include<CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include<CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>
#include<CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_mesh_processing/orientation.h>
#include "lcc_to_polygon_soup.h"
#include "lcc_triangulate_faces.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

using namespace CGAL;
using namespace std;
using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;

namespace PMP=CGAL::Polygon_mesh_processing;

// Définition du LCC 3D
typedef CGAL::Linear_cell_complex_for_combinatorial_map<
    3, 3, CGAL::Linear_cell_complex_traits<3, Kernel>> LCC;

template<typename LCC>
Kernel::Vector_3 triangle_normal(const std::vector<typename LCC::Point>& points,
                                 const std::vector<std::size_t>& triangle)
{
  // Sommets A,B,C
  // Produit vectoriel de AB-> avec AC->
  return CGAL::cross_product(points[triangle[0]] - points[triangle[1]],
                             points[triangle[0]] - points[triangle[2]]);
}

template<typename LCC>
bool is_triangle_external(const std::vector<typename LCC::Point>& points,
                          const std::vector<std::size_t>& triangle,
                          const std::vector<std::vector<std::size_t>>& polygons)
{
  // Barycentre du triangle
  // Deux possibilités : CGAL::barycenter() ou CGAL::centroid() -> différences ?
  typename LCC::Point barycenter = CGAL::centroid(points[triangle[0]], 
                                                  points[triangle[1]], 
                                                  points[triangle[2]]);

  // Normale au triangle
  Kernel::Vector_3 normal = triangle_normal<LCC>(points, triangle);
  
  // Normalisation de la normale
  if(normal.squared_length() > 0)
  {
    normal = normal / std::sqrt(normal.squared_length());
  }

  // Rayon orienté côté positif
  Kernel::Ray_3 pos_ray(barycenter + 1e-6 * normal, normal);
  // Rayon orienté côté négatif
  Kernel::Ray_3 neg_ray(barycenter - 1e-6 * normal, -normal);

  unsigned int pos_ray_intersections = 0;
  unsigned int neg_ray_intersections = 0;

  // Test intersections 
  for(const std::vector<std::size_t>& other : polygons)
  {
    if(other == triangle)
      continue;

    

    if(CGAL::do_intersect(pos_ray, 
                          Kernel::Triangle_3(points[other[0]], points[other[1]], points[other[2]])))
    {
      ++pos_ray_intersections;
    }

    if(CGAL::do_intersect(neg_ray, 
                          Kernel::Triangle_3(points[other[0]], points[other[1]], points[other[2]])))
    {
      ++neg_ray_intersections;
    }
  }

  // Si les deux rayons ont un nombre d'intersections impair -> La face est interne
  if(pos_ray_intersections != 0 && neg_ray_intersections != 0)
    return false;
  else
    return true;
}

// Renvoi le lcc raffiné sans aucune opération booléenne effectuée
void lcc_refinement(const LCC& lcc, LCC& refined_lcc)
{
    LCC lcc_triangulate = lcc;
    // Triangulisation de toutes les 2-cells (faces)
    triangulate_all_faces<LCC>(lcc_triangulate);

    vector<LCC::Point> points;
    vector<vector<size_t>> polygons;

    // Extrait un tableau de 1-cells (points) et un tableau de polygones (faces triangulées)
    lcc_to_polygon_soup<LCC>(lcc_triangulate, points, polygons);

    /* Raffinement géométrique 
     (calcul de toute les intersections entre triangles et produit une nouvelle soupe ou aucun
     triangle n'intersecte l'intérieur d'un autre triangle)
    */
    PMP::autorefine_triangle_soup(points, polygons, parameters::snap_grid_size(32));

    // Reconstruction d'un LCC avec la soupe raffinée
    polygon_soup_to_lcc(points, polygons, refined_lcc, false, true, false);
}

void lcc_union(const LCC& lcc, LCC& union_lcc)
{
    LCC lcc_triangulate = lcc;

    // Triangulisation de toutes les 2-cells (faces)
    triangulate_all_faces<LCC>(lcc_triangulate);

    vector<LCC::Point> points;
    vector<vector<size_t>> polygons;

    // Extrait un tableau de 1-cells (points) et un tableau de polygones (faces triangulées)
    lcc_to_polygon_soup<LCC>(lcc_triangulate, points, polygons);

    /*
    // Test affichage de la soupe 
    cout << "Points (" << points.size() << ") :" << endl;
    for (size_t i = 0; i < points.size(); ++i)
    {
        const auto& p = points[i];
        cout << "  [" << i << "] = (" 
                  << p[0] << ", " << p[1] << ", " << p[2] << ")" << endl;
    }

    cout << "Polygons (" << polygons.size() << ") :" << endl;
    for (size_t i = 0; i < polygons.size(); ++i)
    {
        const auto& poly = polygons[i];
        cout << "  Polygon " << i << " : ";
        for (size_t idx : poly)
            cout << idx << " ";
        cout << endl;
    }

    /* Raffinement géométrique 
     (calcul de toute les intersections entre triangles et produit une nouvelle soupe ou aucun
     triangle n'intersecte l'intérieur d'un autre triangle)
    */
    PMP::autorefine_triangle_soup(points, polygons, parameters::snap_grid_size(32));

    // Supprime les faces intérieures
    //std::vector<std::vector<std::size_t>> triangles_to_erase;
    for(auto it = polygons.begin(); it != polygons.end(); )
    {
    if(!is_triangle_external<LCC>(points, *it, polygons))
    {
        cout << "J'ai trouvé un volume qui n'est pas l'union" << endl;
        it = polygons.erase(it);   // erase retourne l’itérateur suivant
    }
    else
        ++it;
    }

    // Reconstruction d'un LCC avec la soupe raffinée
    polygon_soup_to_lcc(points, polygons, union_lcc, false, true, false);
}

void lcc_intersection(const LCC& lcc, LCC& intersection_lcc)
{
    LCC lcc_triangulate = lcc;

    // Triangulisation de toutes les 2-cells (faces)
    triangulate_all_faces<LCC>(lcc_triangulate);

    vector<LCC::Point> points;
    vector<vector<size_t>> polygons;

    // Extrait un tableau de 1-cells (points) et un tableau de polygones (faces triangulées)
    lcc_to_polygon_soup<LCC>(lcc_triangulate, points, polygons);

    /*
    // Test affichage de la soupe 
    cout << "Points (" << points.size() << ") :" << endl;
    for (size_t i = 0; i < points.size(); ++i)
    {
        const auto& p = points[i];
        cout << "  [" << i << "] = (" 
                  << p[0] << ", " << p[1] << ", " << p[2] << ")" << endl;
    }

    cout << "Polygons (" << polygons.size() << ") :" << endl;
    for (size_t i = 0; i < polygons.size(); ++i)
    {
        const auto& poly = polygons[i];
        cout << "  Polygon " << i << " : ";
        for (size_t idx : poly)
            cout << idx << " ";
        cout << endl;
    }

    /* Raffinement géométrique 
     (calcul de toute les intersections entre triangles et produit une nouvelle soupe ou aucun
     triangle n'intersecte l'intérieur d'un autre triangle)
    */
    PMP::autorefine_triangle_soup(points, polygons, parameters::snap_grid_size(32));

    // Supprime les faces intérieures
    std::vector<decltype(polygons.begin())> to_remove;
    for(auto it = polygons.begin(); it != polygons.end(); )
    {
        if(is_triangle_external<LCC>(points, *it, polygons))
        {
            to_remove.push_back(it);   // erase retourne l’itérateur suivant
        }
        ++it;
    }
    for(auto rit = to_remove.rbegin(); rit != to_remove.rend(); ++rit)
    {
        polygons.erase(*rit);
    }

    // Reconstruction d'un LCC avec la soupe raffinée
    polygon_soup_to_lcc(points, polygons, intersection_lcc, false, true, false);
}

#endif