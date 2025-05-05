#include <iostream>
#include <fstream>
// #include <algorithm>

#include <mfem.hpp>
#include "solver.hpp"  // Solver personnalisé pour le problème


using namespace mfem;
using namespace std;


int main() {

    int vis_steps = 100;
    bool visualization = true;
    int precision = 8;

    real_t t_final = 100000.;
    real_t dt = 0.1;
    real_t t = 0.;

    // real_t rho_ = 1.;
    // real_t c_ = 1.;

    real_t rho_times_c_ = 3.9e6;
    // real_t rho_times_c_ = 1e6;

    // real_t k_ = 0.533;
    real_t k_ = 10;

    unique_ptr<ODESolver> ode_solver = ODESolver::Select(4); // Runge kuta

    int n;       
    int order = 1;  

    const char *path = "../Disque2.msh";
    // const char* path = "D:/Documents/projets/MFEM/mfem/data/periodic-hexagon.mesh";
    Mesh mesh(path, 1, 1);

    mesh.UniformRefinement();
    // mesh.UniformRefinement();

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;


    FiniteElementCollection *fec = new H1_FECollection(order, dim);    
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);


    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;


    GridFunction T(fespace);
    T = 0;

    FunctionCoefficient q_function([&](const Vector &x) { return q_func(x, t); });
    LinearForm q(fespace);
    q.AddDomainIntegrator(new DomainLFIntegrator(q_function));


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
        q_function.SetTime(t);
        // Vector Point(2);
        // Point(0) = 0.4814586436738996;
        // Point(1) = 0.1348983855785123;
        // std::cout << q_function.Eval(Point); 

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
