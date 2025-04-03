#include "mfem.hpp"
#include <iostream>

using namespace std;
using namespace mfem;

int main(){

    const char *path = "../tete.msh";
    Mesh mesh(path);
    mesh.PrintInfo(cout);

    return 0;
}
