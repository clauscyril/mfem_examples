#include <iostream>
#include <fstream>
// #include <algorithm>

#include <mfem.hpp>
#include "solver.hpp"  // Solver personnalisé pour le problème


using namespace mfem;
using namespace std;



int main() {

    int vis_steps = 1;
    bool visualization = true;
    int precision = 8;

    real_t t_final = 3;
    real_t dt = 0.04;
    real_t t = 0.;

    // real_t rho_ = 1.;
    // real_t c_ = 1.;

    real_t rho_times_c_ = 3.9e6;
    // real_t rho_times_c_ = 1e6;

    real_t k_ = 0.533;
    // real_t k_ = 100;

    unique_ptr<ODESolver> ode_solver = ODESolver::Select(4); // Runge kuta
   
    int order = 1;  



    // Mesh Carré généré par mfem
    int nx = 100, ny = 100;
    double sx = 0.032, sy = 0.032;

    Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, true, sx, sy);

    // Centrer le maillage autour de (0,0)
    for (int i = 0; i < mesh.GetNV(); i++) // GetNV() = number of vertices
    {
        double *v = mesh.GetVertex(i);
        v[0] -= sx / 2.0;  // Translation en x
        v[1] -= sy / 2.0;  // Translation en y
    }


    // // Mesh issu de gmsh 
    // const char *path = "../disque.msh";
    // Mesh mesh(path, 1, 1);

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;


    FiniteElementCollection *fec = new H1_FECollection(order, dim);    
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);


    GridFunction q_grid(fespace);


    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;


    compute_Q(mesh,fespace, q_grid);


    GridFunctionCoefficient q_coeff(&q_grid);

    GridFunction T(fespace);
    T = 37;


    LinearForm q(fespace);
    q.AddDomainIntegrator(new DomainLFIntegrator(q_coeff));


    BilinearForm m(fespace);
    BilinearForm k(fespace);

    // ConstantCoefficient alpha(rho_*c_);
    ConstantCoefficient alpha(rho_times_c_);
    m.AddDomainIntegrator(new MassIntegrator(alpha));

    ConstantCoefficient k_coef(-k_);
    k.AddDomainIntegrator(new DiffusionIntegrator(k_coef));

    k.Assemble();
    m.Assemble();

    k.Finalize(0);
    m.Finalize(0);

    q.Assemble(); 

    socketstream sout;

    char vishost[] = "localhost";
    int  visport   = 19916;
    sout.open(vishost, visport);
    if (!sout)
    {
        cout << "Unable to connect to GLVis server at "
            << vishost << ':' << visport << endl;
        visualization = false;
        cout << "GLVis visualization disabled.\n";
    }
    else
    {
        sout.precision(precision);
        sout << "solution\n" << mesh << T << "\nkeys j\n";
        sout << "pause\n";
        sout << flush;
        cout << "GLVis visualization paused."
            << " Press space (in the GLVis window) to resume it.\n";
    }
    



    // Solver
    FE_Evolution adv(m, k, q);

    // real_t t = 0.0;
    adv.SetTime(t);
    ode_solver->Init(adv);



    bool done = false;
    for (int ti = 0; !done; )
    {
        real_t dt_real = min(dt, t_final - t);
        ode_solver->Step(T, t, dt_real);
        ti++;

        done = (t >= t_final - 1e-8*dt);

        adv.SetTime(t);
        // q_coeff.SetTime(t);
        // q.Assemble();
       
        if (done || ti % vis_steps == 0)
        {
            cout << "Time step: " << ti << ", Time: " << t << endl;
        
            if (visualization)
            {
                sout.precision(8);
                sout << "solution\n" << mesh << T
                     << "window_title 'Temperature'"
                     << flush;
                // std::cin.get(); // Décommenter pour pause manuelle
            }
        }
    }

    return 0;
}
