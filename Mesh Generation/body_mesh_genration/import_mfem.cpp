#include "mfem.hpp"
#include <iostream>

using namespace std;
using namespace mfem;

int main() {

    const char *path = "meshs/hand_all.msh";
    Mesh mesh(path, 1, 1);
    
    mesh.PrintInfo(cout);

    return 0;
}

