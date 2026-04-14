#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/draw_linear_cell_complex.h>

using namespace CGAL;
using namespace std;

// Définition du LCC 3D
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;

typedef LCC::Dart_descriptor Dart_descriptor;
typedef LCC::Point           Point;


// TEST MAILLAGES SURFACIQUES
// CREATION D'UN CUBE, METHODE MANUELLE AVEC COMMANDES LCC
// A NOTER QUE LES FACES SONT INDEPENDANTES DONC CA NE DONNE PAS UN VRAI CUBE
int main()
{
    LCC lcc;

    // 2-cell : un carré divisé en deux triangles
    // p3 ---- p2
    //  |       |
    //  |       |
    // p0 ---- p1

    // Fusion de plusieurs 2-cell : 1 cube, 8 points, 12 arrêtes, 6 faces carré (12 triangles)
    
    //      p7 ---- p6
    //    / |      /|
    //   /  |     / |
    //  /   p4 --/- p3 
    // p3 ---- p2  /
    //  |       | /
    //  |       |/
    // p0 ---- p1

    Point p0(0, 0, 0);
    Point p1(1, 0, 0);
    Point p2(1, 1, 0);
    Point p3(0, 1, 0);
    Point p4(0, 0, 1);
    Point p5(1, 0, 1);
    Point p6(1, 1, 1);
    Point p7(0, 1, 1);

    // Face 1 (Avant) : p0, p1, p2, p3
    auto f1_t1 = lcc.make_triangle(p0, p2, p1);
    auto f1_t2 = lcc.make_triangle(p0, p2, p3);
    lcc.sew<2>(f1_t1, f1_t2);

    // Face 2 (Arrière) : p4, p5, p6, p7
    auto f2_t1 = lcc.make_triangle(p4, p6, p5);
    auto f2_t2 = lcc.make_triangle(p4, p6, p7);
    lcc.sew<2>(f2_t1, f2_t2);

    // Face 3 (Gauche) : p0, p3, p4, p7
    auto f3_t1 = lcc.make_triangle(p7, p0, p3);
    auto f3_t2 = lcc.make_triangle(p7, p0, p4);
    lcc.sew<2>(f3_t1, f3_t2);

    // Face 4 (Droite) : p1, p2, p5, p6
    auto f4_t1 = lcc.make_triangle(p6, p1, p2);
    auto f4_t2 = lcc.make_triangle(p6, p1, p5);
    lcc.sew<2>(f4_t1, f4_t2);

    // Face 5 (Haut) : p2, p3, p6, p7
    auto f5_t1 = lcc.make_triangle(p2, p7, p3);
    auto f5_t2 = lcc.make_triangle(p2, p7, p6);
    lcc.sew<2>(f5_t1, f5_t2);

    // Face 6 (Bas) : p0, p1, p4, p5
    auto f6_t1 = lcc.make_triangle(p1, p4, p0);
    auto f6_t2 = lcc.make_triangle(p1, p4, p5);
    lcc.sew<2>(f6_t1, f6_t2);

    
    // AFFICHAGE du LCC avec le viewer CGAL
    draw(lcc);

    return 0;
}