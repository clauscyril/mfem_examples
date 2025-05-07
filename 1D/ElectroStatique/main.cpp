#include "mfem.hpp"
using namespace mfem;
using namespace std;

int main()
{
    // Crée un maillage 1D [0, 1] avec 10 segments
    Mesh mesh = Mesh::MakeCartesian1D(10, 1.0);

    // Attribue un attribut différent aux bords
    // Supposons que le bord gauche (point 0) a l'attribut 1
    // et le bord droit (point 1) a l'attribut 2
    mesh.SetBdrAttribute(0, 1); // point de gauche
    mesh.SetBdrAttribute(1, 2); // point de droite

    // Espace d'éléments finis
    int order = 1;
    H1_FECollection *fec = new H1_FECollection(order, mesh.Dimension());
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);

     // On précise quels noeuds sont des noeuds de Dirichlet
    Array<int> ess_tdof_list;
    for (int i = 0; i < mesh.GetNV(); i++) {
        const double *u = mesh.GetVertex(i);
        std::cout << u[0] << std::endl;
        if (u[0] == 0 or u[0] == 1)
            ess_tdof_list.Append(i);
    }
 
     GridFunction v(fespace);
     v = 0.f;
 
     for (int i = 0; i < ess_tdof_list.Size(); i++){
         const double *u = mesh.GetVertex(ess_tdof_list[i]);
         if (u[0] == 0) 
            v(ess_tdof_list[i]) = 120/sqrt(2);
        else 
            v(ess_tdof_list[i]) = 0;
     }
 

    LinearForm b(fespace);
    // b = 0;
    b.Assemble();

    BilinearForm a(fespace);
    ConstantCoefficient sigma(0.258);
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    
    OperatorPtr A;
    Vector B, V;

    a.FormLinearSystem(ess_tdof_list, v, b, A, V, B);
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, V, 1, 1000, 1e-12, 0.f);

    // On récupère les solutions 
    a.RecoverFEMSolution(B, b, v);

    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock_v(vishost, visport);
    sol_sock_v.precision(8);
    sol_sock_v << "solution\n" << mesh << v
                << "window_title 'Solution: Potentiel'" 
                << "pause\n" << "keys c\n" << std::flush;

    return 0;
}
