#include <iostream>
#include <fstream>
// #include <algorithm>

#include <mfem.hpp>
#include "solver.hpp" // Solver personnalisé pour le problème

using namespace mfem;
using namespace std;

real_t dir_bc_func(const Vector &x);


int main()
{
    int vis_steps = 1;
    bool visualization = true;
    int precision = 8;

    // Paramètres de la simulation en temps
    real_t t_final = 10;
    real_t dt = 0.01;
    real_t t = 0.;

    // Paramètres du système
    real_t rho_times_c_ = 3.9e6;
    real_t k_ = 0.533;

    unique_ptr<ODESolver> ode_solver = ODESolver::Select(4); // Runge kuta

    // Ordre des fonctions 
    int order = 1;

    Mesh *mesh = nullptr;

    const char *path = "../meshs/mesh.msh";
    mesh = new Mesh(path, 1, 1);
   
    int ne = mesh->GetNE();
    int dim = mesh->Dimension();
    int spaceDim = mesh->SpaceDimension();
    
    // real_t z = 8e-3;
    real_t min_x(100);
    real_t min_y(100);
    real_t min_z(100);

    for (int i = 0; i < mesh->GetNV(); i++) // GetNV() = number of vertices
        {
            double *v = mesh->GetVertex(i);
            if (v[0] < min_x)
                min_x = v[0];
            if (v[1] < min_y)
                min_y = v[1];
            if (v[2] < min_z)
                min_z = v[2];
        }

    int id(0);
    int id_1mm(0);

    real_t x = 1.5125;
    real_t y = 3;
    real_t z = 2.25;
    real_t d_min = 10;

    real_t x_1mm = 2.875;
    real_t y_1mm = 3;
    real_t z_1mm = 2.25;
    real_t d_min_1mm = 10;

    for (int i = 0; i < mesh->GetNV(); i++) // GetNV() = number of vertices
    {
        double *v = mesh->GetVertex(i);
        real_t d, d_1mm;
        d = sqrt(pow(v[0] - min_x - x, 2) + pow(v[1] - min_y - y, 2) + pow(v[2] - min_z - z, 2));
        d_1mm = sqrt(pow(v[0] - min_x - x_1mm, 2) + pow(v[1] - min_y - y_1mm, 2) + pow(v[2] - min_z - z_1mm, 2));

        if (d<d_min){
            d_min = d;
            id = i;
        }
        if (d_1mm < d_min_1mm){
            d_min_1mm = d_1mm;
            id_1mm = i;
        }

        v[0] -= min_x;
        v[0] *= 1e-3;  // Translation en x
        v[1] -= min_y;
        v[1] *= 1e-3;  // Translation en y
        v[2] -= min_z;
        v[2] *= 1e-3;
    }

    std::cout << "id milieu : " << id << "\nid 1mm : " << id_1mm << std::endl;

    real_t *u = mesh->GetVertex(id);
    std::cout << "Coordonnées de id : " << u[0] << ", " << u[1] << ", " << u[2] << std::endl;
 
    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;

    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh->bdr_attributes.Size() << endl;

    // Cas 3D reel :
    Array<int> dir_bdr_v(mesh->bdr_attributes.Max());
    dir_bdr_v = 0;
    dir_bdr_v[0] = 1;
    dir_bdr_v[1] = 1;


    Array<int> ess_tdof_list;
    fespace->GetEssentialTrueDofs(dir_bdr_v, ess_tdof_list);
    

    GridFunction v(fespace);
    FunctionCoefficient dir_bc_coef(dir_bc_func);

    std::cout << "test" << std::endl;
    v.ProjectBdrCoefficient(dir_bc_coef, dir_bdr_v);

    // Calcul le potentiel V à partir des conditions aux limites préicsées
    compute_V(*mesh, fespace, ess_tdof_list, v);


    // À partir de ce potentiel, on peut determiner le terme source q
    real_t periode = 0.22 + 0.56;
    real_t t_on = 0.22;
    QCoefficient q_coeff(v, periode, t_on);
    GridFunction q_grid(fespace);
    q_grid.ProjectCoefficient(q_coeff);


    GridFunction T(fespace);
    T = 31;

    real_t T_value = T(id);
    real_t T_1mm = T(id_1mm);

    real_t h = 25;
    real_t T_inf = 25; 
    

    Array<int> robin_bdr(mesh->bdr_attributes.Max());
    robin_bdr = 0; 
    robin_bdr[2] = 1; // Bas 

    Array<int> dirichlet_bdr(mesh->bdr_attributes.Max());
    dirichlet_bdr = 0;
    dirichlet_bdr[3] = 1;

    ConstantCoefficient T_D(31.0);
    T.ProjectBdrCoefficient(T_D, dirichlet_bdr);

    // Marquer les dofs pour exclusion du système
    Array<int> ess_tdof_list_T;
    fespace->GetEssentialTrueDofs(dirichlet_bdr, ess_tdof_list_T);

    ConstantCoefficient h_coeff(-h);
    ConstantCoefficient h_times_T_inf(h*T_inf);

    LinearForm q(fespace);
    q.AddDomainIntegrator(new DomainLFIntegrator(q_coeff));

    q.AddBoundaryIntegrator(new BoundaryLFIntegrator(h_times_T_inf), robin_bdr);

    BilinearForm m(fespace);
    BilinearForm k(fespace);

    // ConstantCoefficient alpha(rho_*c_);
    ConstantCoefficient rho_times_c_coeff(rho_times_c_);
    m.AddDomainIntegrator(new MassIntegrator(rho_times_c_coeff));

    ConstantCoefficient k_coef(-k_);
    k.AddDomainIntegrator(new DiffusionIntegrator(k_coef));
    k.AddBoundaryIntegrator(new BoundaryMassIntegrator(h_coeff), robin_bdr);
    k.Assemble();
    m.Assemble();

    k.Finalize(0);
    m.Finalize(0);

    q.Assemble();

    socketstream sout;
    socketstream sout_q;

    char vishost[] = "localhost";
    int visport = 19916;
    sout.open(vishost, visport);
    sout_q.open(vishost, visport);
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
        sout_q.precision(precision);
        sout << "solution\n" << *mesh << T << "\nkeys j\n" << flush;
        sout_q << "solution\n" << *mesh << q_grid << "\nkeys j\n" << flush;
    }

    // Solver
    FE_Evolution adv(m, k, q);

    adv.SetTime(t);
    ode_solver->Init(adv);

    // Sauvegarde  de la valeur de T entre les électrode dans un fichier :
    std::ofstream data_file("./data.csv");
    data_file << "Time;Temperature;T_1mm\n0;30;30\n";

    bool done = false;
    for (int ti = 0; !done;)
    {
        real_t dt_real = min(dt, t_final - t);
        ode_solver->Step(T, t, dt_real);
        T.ProjectBdrCoefficient(T_D, dirichlet_bdr);
        ti++;

        done = (t >= t_final - 1e-8 * dt);

        adv.SetTime(t);
        q_coeff.SetTime(t);
        q_grid.ProjectCoefficient(q_coeff);
        q.Assemble();

        T_value = T(id);
        T_1mm = T(id_1mm);
        data_file << t << ";" << T_value << ";" << T_1mm << std::endl; 

        // Visualisation avec GLVis
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
                sout_q.precision(precision);
                sout_q << "solution\n" << *mesh << q_grid <<  "window_title 'Terme source q'" << flush;
                // std::cin.get(); // Décommenter pour pause manuelle
            }
        }
    }

    delete mesh;
    return 0;
}


real_t dir_bc_func(const Vector &x){
    real_t V0 = 95/sqrt(2);
    if (x[2] <= 2e-3) {
        return (real_t)0.;
    } else {
        return V0; 
    }
}