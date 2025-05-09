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

    real_t t_final = 10;
    real_t dt = 0.01;
    real_t t = 0.;

    real_t rho_times_c_ = 3.9e6;

    real_t k_ = 0.133;

    real_t V0 = 100/sqrt(2);

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
        const char *path = "../rect3D.msh";
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

    int id(0);
    int id_1mm(0);

    real_t x = 1.5125e-3;
    real_t y = 2.25e-3;
    real_t d_min = 1;

    real_t x_1mm = 2.875e-3;
    real_t y_1mm = 2.25e-3;
    real_t d_min_1mm = 1;

    real_t z = 18e-3;


    // On détermines les conditions aux limites de type Dirichlet du potentiel ici. (Permet de généraliser le fichier solver)
    Array<int> ess_tdof_list;
    for (int i = 0; i < mesh->GetNV(); i++)
    {
        // On cherche d'abord les noeuds les plus proches de là ou l'on souhaite mesurer T
        const double *u = mesh->GetVertex(i);
        real_t d, d_1mm;
        if (dim == 2){
            d = sqrt(pow(u[0] - x,2) + pow(u[1] - y,2));
            d_1mm = sqrt(pow(u[0] - x_1mm,2) + pow(u[1] - y_1mm,2));
        } else {
            d = sqrt(pow(u[0] - x,2) + pow(u[1] - y,2) + pow(u[2] - z, 2));
            d_1mm = sqrt(pow(u[0] - x_1mm,2) + pow(u[1] - y_1mm,2) + pow(u[2] - z, 2));
        }
        if (d < d_min) {
            d_min = d;
            id = i;
        } else if(d_1mm < d_min_1mm) {
            d_min_1mm = d_1mm;
            id_1mm = i;
        }

        // Electrodes carrées
        if (abs(u[0] - 1.5125e-3) < 1.15e-3 and (u[1] <= 1.125e-3 or u[1] >= 3.25e-3))
        {
            if (dim == 3) {
                if (u[2] > 15e-3) {
                    ess_tdof_list.Append(i);
                }
            } else {
                ess_tdof_list.Append(i);
            } 
        }

        // // Électrode circulaire
        // if ( (pow(u[0] - 1.5125e-3, 2) + pow(u[1] - 4.5e-3,2) <=  pow(1.25e-3,2)) or (pow(u[0] - 1.5125e-3, 2) + pow(u[1],2) <=  pow(1.25e-3,2) ))
        // {
        //     ess_tdof_list.Append(i);
        // }
    }

    // Affiche les noeuds les plus proches des positions où l'on souhaite mesurer T
    real_t *u = mesh->GetVertex(id);
    std::cout << "id : " << id << ", x = " << u[0] << ", y = " << u[1] << std::endl;
    u = mesh->GetVertex(id_1mm);
    std::cout << "id : " << id_1mm << ", x = " << u[0] << ", y = " << u[1] << std::endl;

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
    // // 2D
    // robin_bdr[0] = 1; // Bas (attribut 1 → index 0)  
    // robin_bdr[2] = 1; // Haut (attribut 3 → index 2)

    // 3D 
    robin_bdr[3] = 1; // Bas 
    robin_bdr[5] = 1; // Haut

    Array<int> dirichlet_bdr(mesh->bdr_attributes.Max());
    dirichlet_bdr = 0;
    dirichlet_bdr[4] = 1;

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

    // real_t t = 0.0;
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
