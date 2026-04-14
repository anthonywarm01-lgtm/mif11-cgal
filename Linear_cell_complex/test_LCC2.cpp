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
// CREATION D'UN CUBE, METHODE AUTOMATIQUE AVEC LE MAKE HEXAHEDRON de LCC
int main()
{

    // Fusion de 1-cell : 1 cube, 8 points, 12 arrêtes, 6 faces carré (12 triangles ?)
    
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

    // Création du LCC
    LCC lcc;

    lcc.make_hexahedron(p0, p1, p2, p3, p7, p4, p5, p6);

    // Affichage rapide pour vérification, doit afficher 24 darts (12 arrêtes * 2)
    lcc.display_characteristics(cout) << ", valid=" << lcc.is_valid() << endl;

    // AFFICHAGE du LCC avec le viewer CGAL
    draw(lcc);

    return 0;
}