#include "mfem.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
    const char *mesh_file = "D:/Documents/projets/MFEM/mfem_examples/2D/serial/data/test.msh";
    ifstream mesh_stream(mesh_file);

    if (!mesh_stream)
    {
        cerr << "Error: Cannot open mesh file!" << endl;
        return 1;
    }

    // Try reading the mesh
    try
    {
        Mesh mesh(mesh_file, 1, 1);
        cout << "Mesh loaded successfully!" << endl;
        mesh.PrintInfo(cout);
    }
    catch (exception &e)
    {
        cerr << "Error reading mesh: " << e.what() << endl;
    }

    return 0;
}
