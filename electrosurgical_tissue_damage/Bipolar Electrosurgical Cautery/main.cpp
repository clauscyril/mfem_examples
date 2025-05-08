#include <iostream>
#include <fstream>
// #include <algorithm>

#include <mfem.hpp>
#include "solver.hpp" // Solver personnalisé pour le problème

using namespace mfem;
using namespace std;

int main()
{

    int vis_steps = 1;
    bool visualization = true;
    int precision = 8;

    real_t t_final = 3;
    real_t dt = 0.01;
    real_t t = 0.;

    real_t rho_times_c_ = 3.9e6;

    real_t k_ = 0.533;

    real_t V0 = 75;
    real_t alpha = 0.531;

    unique_ptr<ODESolver> ode_solver = ODESolver::Select(4); // Runge kuta

    int order = 1;

    bool generateMesh = 0;

    Mesh *mesh = nullptr;

    real_t lc = 4.5e-3 / 50;
    double sx = 10e-3, sy = 4.5e-3;

    if (generateMesh)
    {
        // MESH 2D
        // Mesh Carré généré par mfem
        int nx = floor(sx / lc), ny = floor(sy / lc);
        mesh = new Mesh(Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, true, sx, sy));
    }
    else
    {
        // Mesh issu de gmsh
        const char *path = "../rect.msh";
        mesh = new Mesh(path, 1, 1);
    }

    int ne = mesh->GetNE();
    int dim = mesh->Dimension();
    int spaceDim = mesh->SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;

    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh->bdr_attributes.Size() << endl;

    int id, id_1mm;

    // On détermines les conditions aux limites de type Dirichlet du potentiel ici. (Permet de généraliser le fichier solver)
    Array<int> ess_tdof_list;
    for (int i = 0; i < mesh->GetNV(); i++)
    {
        const double *u = mesh->GetVertex(i);
        // std::cout << i << " : " << u[0] << ", " << u[1] << std::endl;

        // Electrodes carrées
        if (abs(u[0] - 1.5125e-3) < 1.15e-3 and (u[1] <= 1.125e-3 or u[1] >= 3.25e-3))
        {
            ess_tdof_list.Append(i);
        }

        // // Électrode circulaire
        // if ( (pow(u[0] - 1.5125e-3, 2) + pow(u[1] - 4.5e-3,2) <=  pow(1.25e-3,2)) or (pow(u[0] - 1.5125e-3, 2) + pow(u[1],2) <=  pow(1.25e-3,2) ))
        // {
        //     ess_tdof_list.Append(i);
        // }
    }

    GridFunction v(fespace);
    v = 0.f;

    for (int i = 0; i < ess_tdof_list.Size(); i++)
    {
        const double *u = mesh->GetVertex(ess_tdof_list[i]);
        if (u[1] > sy / 2)
        {
            v(ess_tdof_list[i]) = V0;
        }
        else
        {
            v(ess_tdof_list[i]) = 0.f;
        }
    }

    // GridFunction q_grid(fespace);

    // Calcul le potentiel V à partir des conditions aux limites préicsées
    compute_V(*mesh, fespace, ess_tdof_list, v);

    // À partir de ce potentiel, on peut determiner le terme source q
    real_t periode = 0.22 + 0.56;
    real_t t_on = 0.22;
    QCoefficient q_coeff(v, periode, t_on);

    GridFunction T(fespace);
    T = 30;

    LinearForm q(fespace);
    q.AddDomainIntegrator(new DomainLFIntegrator(q_coeff));

    BilinearForm m(fespace);
    BilinearForm k(fespace);

    // ConstantCoefficient alpha(rho_*c_);
    ConstantCoefficient rho_times_c_coeff(rho_times_c_);
    m.AddDomainIntegrator(new MassIntegrator(rho_times_c_coeff));

    ConstantCoefficient k_coef(-k_);
    k.AddDomainIntegrator(new DiffusionIntegrator(k_coef));

    k.Assemble();
    m.Assemble();

    k.Finalize(0);
    m.Finalize(0);

    q.Assemble();

    socketstream sout;

    char vishost[] = "localhost";
    int visport = 19916;
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
        sout << "solution\n"
             << *mesh << T << "\nkeys j\n";
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
    for (int ti = 0; !done;)
    {
        real_t dt_real = min(dt, t_final - t);
        ode_solver->Step(T, t, dt_real);
        ti++;

        done = (t >= t_final - 1e-8 * dt);

        adv.SetTime(t);
        q_coeff.SetTime(t);
        q.Assemble();

        if (done || ti % vis_steps == 0)
        {
            cout << "Time step: " << ti << ", Time: " << t << endl;

            if (visualization)
            {
                sout.precision(8);
                sout << "solution\n"
                     << *mesh << T
                     << "window_title 'Temperature'"
                     << flush;
                // std::cin.get(); // Décommenter pour pause manuelle
            }
        }
    }

    delete mesh;
    return 0;
}
