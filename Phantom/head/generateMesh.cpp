#include <CGAL/Simple_cartesian.h>
// #include "mfem.hpp"
#include <iostream>

using namespace std;
// using namespace mfem;

int main(){

    // const char *path = "../tete.msh";
    // Mesh mesh(path);
    // mesh.PrintInfo(cout);

    typedef CGAL::Simple_cartesian<double> Kernel;
    Kernel::Point_2 p(1, 2), q(3, 4);
    cout << "Distance: " << CGAL::squared_distance(p, q) << endl;


    return 0;
}
