#include <mfem.hpp>
#include <iostream>

using namespace mfem;
using namespace std;


int main(){

    const char *mesh_file = "../beam-tet.mesh";
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();
    int sdim = mesh->SpaceDimension();

    cout << "Dimension : " << dim << "\nSpace dimension : " << sdim << endl;
    
    int test = 17;
    int *ptr = &test;

    const char *str = "Hello";

    cout << sizeof(ptr) << endl;
    cout << M_PI << endl;
    cout << sizeof(str) << endl;

    // cout << i;
    // cout << static_cast<const void*>(str+1) << endl;   
    return 0;
}