#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include "lcc_boolean_operations.h"

using namespace CGAL;
using namespace std;

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;

typedef LCC::Dart_descriptor Dart_descriptor;
typedef LCC::Point           Point;

// TEST OPERATIONS BOOLEENNES
// Arguments de l'appel : ../../../../../Data/data/meshes/forme_voulu.off
int main()
{
    LCC lcc;
    if(!load_off(lcc, data_file_path("meshes/tetrahedron.off").c_str()))
    {
        cout << "ERROR reading file tetrahedron.off" << endl;
        exit(EXIT_FAILURE);
    }

    if(!load_off(lcc, data_file_path("meshes/sphere.off").c_str()))
    {
        cout << "ERROR reading file sphere.off" << endl;
        exit(EXIT_FAILURE);
    }

    lcc.display_characteristics(std::cout) << ", valid=" << lcc.is_valid() << std::endl;
    draw(lcc);

    /*LCC refined_lcc;
    lcc_refinement(lcc, refined_lcc);

    refined_lcc.display_characteristics(std::cout) << ", valid=" << refined_lcc.is_valid() << std::endl;
    draw(refined_lcc);*/

    LCC union_lcc;
    lcc_union(lcc, union_lcc);

    union_lcc.display_characteristics(std::cout) << ", valid=" << union_lcc.is_valid() << std::endl;
    draw(union_lcc);

    return 0;
}