#include "solver_TD_FOM.hpp"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"
#include <iostream>
#include <cmath>
#include <unordered_set>

using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

void TD_sim_offline(Mesh &mesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization, bool save_snapshot){

    // ****************** Parameters definitions *********************
    real_t Ts = t_f/(num_steps);

    // geometric parameters
    real_t Ri = 9.6e-3/2.0;
    real_t height = 7.59e-3;
    real_t w = 5.3e-3;
    real_t Rout = Ri + w;
    real_t Rm = (Rout - Ri) / 2;
    int iter = 0;

    // ****** Ferrite parameters ******
    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;
    real_t mu = ferrite.mu; 
    
    real_t tau = 1/(2 * M_PI * 1.8e6);

    // int num_steps = 1000;
    real_t t = 0.0;

    real_t A1 = (rho * (sigma * Ts + 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A2 = (rho * (sigma * Ts - 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A3 = -(Ts*sigma -2*eps) / (2*eps + Ts*sigma);

    real_t B1 = 2*mu / (Ts + 2*tau);
    real_t B2 = -(2*mu) / (Ts + 2*tau);
    real_t B3 = -(Ts - 2*tau) / (Ts + 2*tau);

    real_t C1 = Ts*mu/(Ts+2*tau);
    real_t C2 = Ts*mu/(Ts+2*tau);
    real_t C3 = -(Ts-2*tau)/(Ts+2*tau);

    real_t NI = 0; // For now, initialized to 0;


    int order = 1;                      // elements order : 2
    int dim = mesh.Dimension();        // Mesh dimension : 2
    int id = 0;

    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();


    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // Boundary conditions function
    auto bdr_func = [&](const Vector &x, real_t t) -> real_t
    {
        real_t r = x(0); // Accès au premier élément du vecteur x
        return NI_func(t) / (2 * M_PI * r);
    };

    // // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];
    };


    // Function for the correcting factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };
    

    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(&pmesh, fec);
    ParFiniteElementSpace *fespace_E = new ParFiniteElementSpace(&pmesh, fec, dim);

    HYPRE_Int size = fespace->GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of finite element unknowns: " << size << endl;
    }

    ParGridFunction Hn(fespace);     // H^n
    ParGridFunction Hnm1(fespace);   // H^{n-1}
    ParGridFunction db_dt(fespace);
    ParGridFunction db_dtm1(fespace);

    ParGridFunction Bn(fespace);
    ParGridFunction Bnm1(fespace);

    ParGridFunction En(fespace_E);
    ParGridFunction Enm1(fespace_E);

    ParGridFunction Jn(fespace_E);
    ParGridFunction Jn_x(fespace);  
    ParGridFunction Jn_y(fespace);  
    ParGridFunction Jnm1(fespace_E);
    // CurlCustomCoefficient *curl_H_coeff = new CurlCustomCoefficient(&Hn);

    ParGridFunction Pn(fespace);

    // Initializing all values to 0
    Hn = 0;
    Bn = 0;
    Bnm1 = 0;
    En = 0;
    Jn = 0;
    Jnm1 = 0;
    Enm1 = 0;
    Enm1 = 0;

    Hnm1 = 0;  
    db_dt = 0;


    // Definition of the Boundary conditions  :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(pmesh.bdr_attributes.Max());
    dir_bdr = 1; // All the borders have boundary conditions
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);

    
    FunctionCoefficient bdr_coeff([&](const Vector &x) { return bdr_func(x, t); });
    Hnm1.ProjectBdrCoefficient(bdr_coeff, dir_bdr);

    // **** Coefficients for bilinear forms ****
    FunctionCoefficient r_coeff(r_coeff_func);
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);

    ConstantCoefficient A1_coeff(A1);
    ProductCoefficient A1_r(A1_coeff, r_coeff);

    ConstantCoefficient B1_coeff(B1);
    ProductCoefficient B1_r(B1_coeff, r_coeff);

    ProductCoefficient A1_r_inv(A1_r, inv_r_square_coeff);


    ConstantCoefficient A2_coeff(-A2);
    ProductCoefficient A2_r(A2_coeff, r_coeff);

    ConstantCoefficient B2_coeff(-B2);
    ProductCoefficient B2_r(B2_coeff, r_coeff);

    ProductCoefficient A2_r_inv(A2_r, inv_r_square_coeff);

    ConstantCoefficient k(-(B3-A3));
    ProductCoefficient k_r(k, r_coeff);
    // **************************************************

    // ******* Definition of Bilinear and Linear Forms ********
    ParBilinearForm r1(fespace); 
    r1.AddDomainIntegrator(new DiffusionIntegrator(A1_r));
    r1.AddDomainIntegrator(new MassIntegrator(B1_r));
    r1.AddDomainIntegrator(new MassIntegrator(A1_r_inv));
    r1.Assemble();

    ParBilinearForm r2(fespace); 
    r2.AddDomainIntegrator(new DiffusionIntegrator(A2_r));
    r2.AddDomainIntegrator(new MassIntegrator(B2_r));
    r2.AddDomainIntegrator(new MassIntegrator(A2_r_inv));
    r2.Assemble();

    ParBilinearForm m(fespace);
    m.AddDomainIntegrator(new MassIntegrator(k_r));
    m.Assemble();

    ParLinearForm lf_surface(fespace);
    ConstantCoefficient one_coeff(1.);
    lf_surface.AddDomainIntegrator(new DomainLFIntegrator(one_coeff));
    lf_surface.Assemble();

    r1.Finalize();
    r2.Finalize();
    m.Finalize();
    

    // Parameters for the linear system
    // SparseMatrix &A = r1.SpMat();
    Vector rhs(fespace->GetTrueVSize());
    Vector R2Hnm1(fespace->GetTrueVSize());
    Vector MdBdt(fespace->GetTrueVSize());

    Vector x(fespace->GetTrueVSize()); // For defining the linear system, considering Dirichlet's bdc

    // Solver for the linear system
    // UMFPackSolver solver;
    // solver.SetOperator(A);
    CGSolver solver(MPI_COMM_WORLD);
    solver.SetRelTol(1e-12);
    solver.SetMaxIter(2000);
    solver.SetPrintLevel(0);




    // Sockets stream for showing the results using glvis (if visualization is set to true) 
    socketstream sout;
    socketstream sout_e;
    socketstream sout_j;
    socketstream sout_b;

    char vishost[] = "localhost";
    int visport = 19916;

    if (visualization)  // Initialize the ploting of the results
    {
        sout.open(vishost, visport);
        sout_e.open(vishost, visport);
        sout_j.open(vishost, visport);
        sout_b.open(vishost, visport);

        sout << "parallel " << num_procs << " " << myid << "\n";
        sout.precision(8);
        sout << "solution\n" << pmesh << Hn << "\nkeys j\n" << std::flush;
        sout_e << "parallel " << num_procs << " " << myid << "\n";
        sout_e.precision(8);
        sout_e << "solution\n" << pmesh << En << "\nkeys j\n" << std::flush;
        sout_j.precision(8);
        sout_j << "solution\n" << pmesh << Jn << "\nkeys j\n" << std::flush;
        sout_b.precision(8);
        sout_b << "solution\n" << pmesh << Bn << "\nkeys j\n" << std::flush;
    }

    // File for saving the values of the power and flux for each iterations
    std::string name = "./data/fom/TD_" + std::to_string(myid) + ".csv";   
    std::ofstream data_file(name);                         
    data_file << "t;p_eddy;flux;NI\n0;0;0;0\n";  // Intialising the file with coluns names and first values to 0


    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "data/basis";
    const std::string basisFileName = basisName + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasis;
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    int numRowRB, numColumnRB;
    StopWatch solveTimer, assembleTimer, mergeTimer;

    const std::string basisNameEx = "data/basisEx";
    const std::string basisNameEy = "data/basisEy";
    const std::string basisFileNameEx = basisNameEx + std::to_string(id);
    const std::string basisFileNameEy = basisNameEy + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasisEx;
    std::unique_ptr<const CAROM::Matrix> spatialbasisEy;
    CAROM::BasisGenerator *generatorEx;
    CAROM::BasisGenerator *generatorEy;

    // options = new CAROM::Options(newSize, max_num_snapshots, update_right_SV);
    options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots, update_right_SV);

    generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);

    generatorEx = new CAROM::BasisGenerator(*options, isIncremental, basisFileNameEx);
    generatorEy = new CAROM::BasisGenerator(*options, isIncremental, basisFileNameEy);


    int vis_steps = 10;  // Number of iterations between two visualizations
    // ***************  Time iterations *****************
    for (int step = 0; step < num_steps; step++)
    {   
        iter++;

        // std::cout << iter%100 << std::endl;
        t += Ts;                    // New time 
        bdr_coeff.SetTime(t);       // Update boundary conditions fonction's time parameters
        Hn.ProjectBdrCoefficient(bdr_coeff, dir_bdr); // Setting new boundary conditions on H

        // Mult M * dBdt → MdBdt
        m.Mult(db_dtm1, MdBdt);

        // Mult R2 * Hnm1 → R2Hnm1
        r2.Mult(Hnm1, R2Hnm1);

        // RHS = R2Hnm1 + MdBdt = R2* * Hnm1 + M * dBnm1/dt
        rhs = R2Hnm1;
        rhs += MdBdt;
        
        HypreParMatrix A_sys;
        Vector X, B;

        // Form linear system taking boundary conditions into account
        r1.FormLinearSystem(ess_tdof_list, Hn, rhs, A_sys, X, B);

        // Solving Linear system : A_sys * X = B
        solver.SetOperator(A_sys);
        solver.Mult(B, X);
        // Get solution into Hn
        r1.RecoverFEMSolution(X, rhs, Hn);


        // *** Update of GridFunctions from solution Hn ***

        // Update of db_dt with : (db_dt)_n = B1*Hn + B2*Hn-1 + B3*(db_dt)_n-1
        db_dtm1 *= B3;

        ParGridFunction H_temp(fespace); // temporary GridFunction to avoid modifying Hn
        H_temp = Hn;
        H_temp *= B1;  // Hn*B1

        Hnm1 *= B2;     // Hn-1 * B2

        db_dt = 0;
        db_dt += H_temp;
        db_dt += Hnm1;
        db_dt += db_dtm1;  // db_dt up to date

        Hnm1 /= B2;         // Restoring Hnm1 and db_dtm1 for the update of Bn
        db_dtm1 /= B3;
        H_temp /= B1;
    
        // ***** Update of Bn ******
        // Bn = C1 Hn + C2 Hnm1 + C3 Bnm1
        H_temp *= C1;
        Hnm1 *= C2;
        Bnm1 *= C3;

        Bn = H_temp;
        Bn += Hnm1;
        Bn += Bnm1; // Bn up to date

        // Update of n-1 iterations 
        Hnm1 = Hn;          // Storing Hn into Hn-1 for next iteration
        db_dtm1 = db_dt;    // Same for db_dt
        Bnm1 = Bn;          // Same for Bn

        // ***** computation of magnetic flux ***********
        GridFunctionCoefficient B_coeff(&Bn);
        ParGridFunction one(fespace);
        one = 1;
        ParLinearForm lf_flux(fespace);
        lf_flux.AddDomainIntegrator(new DomainLFIntegrator(B_coeff));
        lf_flux.Assemble();

        real_t flux = lf_flux(one);

        // Computation of J and E
        CurlCustomCoefficient J_coeff(&Hn);
        Jn.ProjectCoefficient(J_coeff);  // Jn = curl(Hn)

        ParGridFunction Jn_temp(fespace_E);
        Jn_temp = Jn;
        Jn_temp *=A1;
        Jnm1 *= A2;
        Enm1 *= A3;
        En = 0;
        En+=Jn_temp;
        En+=Jnm1;
        En+=Enm1;

        Enm1 = En;
        Jnm1 = Jn;

        VectorGridFunctionCoefficient E_coeff(&En);
        PowerLossCoefficient_TD P_eddy_coeff(fespace_E, J_coeff, E_coeff);
        ParLinearForm lf(fespace);
        lf.AddDomainIntegrator(new DomainLFIntegrator(P_eddy_coeff));
        lf.Assemble();
        real_t P_eddy_tot = 2*M_PI*lf(one);
        real_t vol = M_PI * height * (Rout*Rout - Ri*Ri);
        real_t P_eddy = P_eddy_tot/vol;

        Pn.ProjectCoefficient(P_eddy_coeff);


        data_file << t << ";" << P_eddy << ";" << flux << ";" << NI_func(t) << std::endl;

        // Visualisation avec GLVis
        if (step % vis_steps == 0 )
        {   
            std::cout << "Time step: " << step << ", Time: " << t << std::endl;


            std::string name2 = "./data/fom/Hn_" + std::to_string(step) + std::to_string(myid)+ ".csv";   
            std::cout << name2 << std::endl;
            std::ofstream data_file2(name2);                         
            data_file2 << "Hn\n";  // Intialising the file with coluns names and first values to 0
            for (int i = 0; i<Hn.Size(); i++) {
                data_file2 << Hn(i) << std::endl;
            }

            if (save_snapshot)
            {
                for (int i = 0; i<Hn.Size(); i++) {
                    Jn_x(i) = Jn(i);
                    Jn_y(i) = Jn(i+Hn.Size());
                }
                
                bool addSample = generator->takeSample(Hn.GetData());
                bool addSampleEx = generatorEx->takeSample(Jn_x.GetData());
                bool addSampleEy = generatorEy->takeSample(Jn_y.GetData());
            }


            if (visualization)
            {
                sout.precision(8);
                sout << "solution\n"
                     << pmesh << Jn_y
                     << "window_title 'Champ H'"
                     << std::flush;
                sout_e << "solution\n" << pmesh << En <<  "window_title 'Champ E'" << std::flush;
                sout_j << "solution\n" << pmesh << Jn <<  "window_title 'Courant J'" << std::flush;
                sout_b << "solution\n" << pmesh << Jn_x <<  "window_title 'Field B'" << std::flush;

                // std::cin.get(); // Décommenter pour pause manuelle
            }
        }
        
    }
    if (save_snapshot) {
        generator->writeSnapshot();
        generatorEx->writeSnapshot();
        generatorEy->writeSnapshot();
        generator->endSamples();
        generatorEx->endSamples();
        generatorEy->endSamples();
    }
    delete generator;
    delete options;
    data_file.close();
    Mpi::Finalize();

}

void TD_sim_online(Mesh &mesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization){
    
    real_t Ts = t_f/(num_steps);

    // geometric parameters
    real_t Ri = 9.6e-3/2.0;
    real_t height = 7.59e-3;
    real_t w = 5.3e-3;
    real_t Rout = Ri + w;
    real_t Rm = (Rout - Ri) / 2;
    int iter = 0;

    // ****** Ferrite parameters ******
    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;
    real_t mu = ferrite.mu; 
    
    real_t tau = 1/(2 * M_PI * 1.8e6);

    // int num_steps = 1000;
    real_t t = 0.0;

    real_t A1 = (rho * (sigma * Ts + 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A2 = (rho * (sigma * Ts - 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A3 = -(Ts*sigma -2*eps) / (2*eps + Ts*sigma);

    real_t B1 = 2*mu / (Ts + 2*tau);
    real_t B2 = -(2*mu) / (Ts + 2*tau);
    real_t B3 = -(Ts - 2*tau) / (Ts + 2*tau);

    real_t C1 = Ts*mu/(Ts+2*tau);
    real_t C2 = Ts*mu/(Ts+2*tau);
    real_t C3 = -(Ts-2*tau)/(Ts+2*tau);

    real_t NI = 0; // For now, initialized to 0;


    int order = 1;                      // elements order : 2
    int dim = mesh.Dimension();        // Mesh dimension : 2
    int id = 0;

    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();


    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // Boundary conditions function
    auto bdr_func = [&](const Vector &x, real_t t) -> real_t
    {
        real_t r = x(0); // Accès au premier élément du vecteur x
        return NI_func(t) / (2 * M_PI * r);
    };

    // // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];
    };


    // Function for the correcting factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };
    


    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(&pmesh, fec);
    ParFiniteElementSpace *fespace_E = new ParFiniteElementSpace(&pmesh, fec, dim);

    HYPRE_Int size = fespace->GlobalTrueVSize();

    ParGridFunction Hn(fespace);     // H^n
    ParGridFunction Hnm1(fespace);   // H^{n-1}
    ParGridFunction Pn(fespace);     // Pn
    ParGridFunction Jnx(fespace);     // Pn
    ParGridFunction Jny(fespace);     // Pn
    // ParGridFunction Jnx_test(fespace);     // Pn

    ParGridFunction Bn(fespace);

    ParGridFunction En(fespace_E);
    ParGridFunction Enm1(fespace_E);

    ParGridFunction Jn(fespace_E);
    ParGridFunction Jnm1(fespace_E);

    // Initializing all values to 0
    Hn = 0;
    Hnm1 = 0;  
    Pn = 0;

    Bn = 0;
    En = 0;
    Enm1 = 0;
    Jn = 0;
    Jnm1 = 0;



    // Definition of the Boundary conditions  :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(pmesh.bdr_attributes.Max());
    dir_bdr = 1; // All the borders have boundary conditions
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);


    FunctionCoefficient bdr_coeff([&](const Vector &x) { return bdr_func(x, t); });
    Hnm1.ProjectBdrCoefficient(bdr_coeff, dir_bdr);  // Not necessary as the initial bdc is 0

    // **** Coefficients for bilinear forms ****
    FunctionCoefficient r_coeff(r_coeff_func);
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);

    ConstantCoefficient A1_coeff(A1);
    ProductCoefficient A1_r(A1_coeff, r_coeff);

    ConstantCoefficient B1_coeff(B1);
    ProductCoefficient B1_r(B1_coeff, r_coeff);

    ProductCoefficient A1_r_inv(A1_r, inv_r_square_coeff);


    ConstantCoefficient A2_coeff(-A2);
    ProductCoefficient A2_r(A2_coeff, r_coeff);

    ConstantCoefficient B2_coeff(-B2);
    ProductCoefficient B2_r(B2_coeff, r_coeff);

    ProductCoefficient A2_r_inv(A2_r, inv_r_square_coeff);

    ConstantCoefficient k(-(B3-A3));
    ProductCoefficient k_r(k, r_coeff);
    // **************************************************


    ParBilinearForm r1(fespace); 
    r1.AddDomainIntegrator(new DiffusionIntegrator(A1_r));
    r1.AddDomainIntegrator(new MassIntegrator(B1_r));
    r1.AddDomainIntegrator(new MassIntegrator(A1_r_inv));
    r1.Assemble();

    ParBilinearForm r2(fespace); 
    r2.AddDomainIntegrator(new DiffusionIntegrator(A2_r));
    r2.AddDomainIntegrator(new MassIntegrator(B2_r));
    r2.AddDomainIntegrator(new MassIntegrator(A2_r_inv));
    r2.Assemble();

    ParBilinearForm m(fespace);
    m.AddDomainIntegrator(new MassIntegrator(k_r));
    m.Assemble();



    r1.Finalize();
    r2.Finalize();
    m.Finalize();
    // curl_x.Finalize();
    // curl_y.Finalize();

    // Parameters for the linear system
    // SparseMatrix &A = r1.SpMat();
    // // HypreParMatrix &A = r1.Ma
    // SparseMatrix &R2 = r2.SpMat();
    // SparseMatrix &M = m.SpMat();
    HypreParMatrix A, R2, M;
    ParLinearForm b(fespace);
    b = 0;
    b.Assemble();
    Vector X, B;
    Vector zero_vector(Hn.Size());
    
    Vector rhs(fespace->GetTrueVSize());
    Vector R2Hnm1(fespace->GetTrueVSize());
    Vector MdBdt(fespace->GetTrueVSize());


    Array<int> ess_tdof_list_fake; // Used to generate the Matrixes of the bilinear forms r2 and m 

    r1.FormLinearSystem(ess_tdof_list, Hn, b, A, X, B); // Here using FormLinearSystem to initialize X and B 
    r2.FormSystemMatrix(ess_tdof_list_fake, R2); // Generating the matrixes R2 and M from the bilinear forms
    m.FormSystemMatrix(ess_tdof_list_fake, M);   // --------------------------------



    // ---------- Generation of curl matrixes -----------------
    std::cout << "test" << std::endl;
    HypreParMatrix *test;
    ParDiscreteLinearOperator curl(fespace, fespace_E);
    CurlCustomCoefficient grad_phi_i(&Hn);
    curl.AddDomainInterpolator(new VectorScalarProductInterpolator(grad_phi_i));

    curl.Assemble();
    curl.Finalize();
    curl.PrintMatlab(std::cout);
    std::cout << curl.Height() << ", " << curl.Width() << std::endl;

    // SparseMatrix Cx(Hn.Size(), Hn.Size()); SparseMatrix Cy(Hn.Size(), Hn.Size());
    // GridFunction grad_phi_x(fespace);          Vector grad_phi_y(Hn.Size());
    // GridFunction curl_phi(fespace_E);

    // for (int i = 0; i<Hn.Size(); i++) {
    //     Pn = 0;
    //     Pn(i) = 1;
    //     CurlCustomCoefficient grad_phi_i(&Pn);
    //     // GradientGridFunctionCoefficient grad_phi_i(&Pn);
    //     curl_phi.ProjectCoefficient(grad_phi_i);
    //     std::cout << "Constructing curl matrixes : " << 100*(real_t)i/Hn.Size() << "%" << std::endl;
    //     for (int j = 0; j<Hn.Size(); j++) {
            
    //         grad_phi_x(i) = curl_phi(i);
    //         if (curl_phi(j) != 0) {
    //             Cx.Add(j, i, curl_phi(j));
    //         }
    //         if (curl_phi(j+Hn.Size()) != 0) {
    //             Cy.Add(j, i, curl_phi(j+Hn.Size()));
    //         }
    //     }
    // }
    // std::cout << "test" << std::endl;
    // ----------------------------------------------------------


    
    // File for saving the values of the power and flux for each iterations
    std::string name = "./data/reduced/TD_" + std::to_string(1) + ".csv";   
    std::ofstream data_file(name);                         
    data_file << "t;p_eddy;flux;NI;flux_reduced\n0;0;0;0;0\n";  // Intialising the file with coluns names and first values to 0



    // Elements for generating the reduced system (LibRom elements)
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "data/basis";
    const std::string basisFileName = basisName + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasis;
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    int numRowRB, numColumnRB;
    StopWatch solveTimer, assembleTimer, mergeTimer;

    const std::string basisNameEx = "data/basisEx";
    const std::string basisNameEy = "data/basisEy";
    const std::string basisFileNameEx = basisNameEx + std::to_string(id);
    const std::string basisFileNameEy = basisNameEy + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasisEx;
    std::unique_ptr<const CAROM::Matrix> spatialbasisEy;
    CAROM::BasisGenerator *generatorEx;
    CAROM::BasisGenerator *generatorEy;



    assembleTimer.Start();  // Timer for estimating the assembling time of the reduced system
    CAROM::BasisReader reader(basisFileName);  // Get Snapshot in a BasisReader object
    CAROM::BasisReader readerEx(basisFileNameEx);  // Get Snapshot in a BasisReader object
    CAROM::BasisReader readerEy(basisFileNameEy);  // Get Snapshot in a BasisReader object
    spatialbasis = reader.getSpatialBasis(); // Generate the reduced basis
    spatialbasisEx = readerEx.getSpatialBasis(); // Generate the reduced basis
    spatialbasisEy = readerEy.getSpatialBasis(); // Generate the reduced basis
    numRowRB = spatialbasis->numRows();
    numColumnRB = spatialbasis->numColumns();
    if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                numColumnRB);

    // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
    DenseMatrix *reducedBasisT = new DenseMatrix(spatialbasis->getData(),
            numColumnRB, numRowRB);
    
    
    // Save the spatial basis in a file
    spatialbasis->print("data/phi");
    spatialbasisEx->print("data/phiEx");
    spatialbasisEy->print("data/phiEy");
    // 21. form ROM operators

    CAROM::Matrix ReducedA(numColumnRB, numColumnRB, false);
    CAROM::Matrix invReducedA(numColumnRB, numColumnRB, false);
    CAROM::Matrix ReducedR2(numColumnRB, numColumnRB, false);
    CAROM::Matrix ReducedM(numColumnRB, numColumnRB, false);
    // CAROM::Matrix ReducedCx(numColumnRB, numColumnRB, false);
    // CAROM::Matrix ReducedCy(numColumnRB, numColumnRB, false);


    // CAROM::Matrix Cx_Carom(Hn.Size(), Hn.Size(), false);
    // CAROM::Matrix Cy_Carom(Hn.Size(), Hn.Size(), false);



    ComputeCtAB(A, *spatialbasis, *spatialbasis, ReducedA);    // Reduced A = phi^t * A * phi
    ComputeCtAB(A, *spatialbasis, *spatialbasis, invReducedA); // inverse of Reduced A 
    invReducedA.inverse();
    ComputeCtAB(R2, *spatialbasis, *spatialbasis, ReducedR2); // Reduced R2
    ComputeCtAB(M, *spatialbasis, *spatialbasis, ReducedM);   // Reduced M
    invReducedA.print("data/A_inv");
    ReducedA.print("data/A");

    // ComputeCurl(Mass, Kx, Cx_Carom);
    // ComputeCurl(Mass, Ky, Cy_Carom);

    // Cy_Carom.print("data/reduced/test_Cy");
    // ComputeCtAB(Cx_Carom, *spatialbasis, *spatialbasis, ReducedA);    // Reduced A = phi^t * A * phi
    // ComputeCtAB(A, *spatialbasis, *spatialbasis, invReducedA); // inverse of Reduced A 

    CAROM::Vector B_carom(B.GetData(), B.Size(), true, false);  // B Compatible with libRom
    CAROM::Vector X_carom(X.GetData(), X.Size(), true, false);  // X Compatible with libRom 
    CAROM::Vector Bn_carom(Bn.GetData(), Bn.Size(), true, false);  // X Compatible with libRom
    
    // CAROM::Vector Ex_carom(zero_vector.GetData(), zero_vector.Size(), true, false);  //
    // CAROM::Vector Ey_carom(zero_vector.GetData(), zero_vector.Size(), true, false);

    // CAROM::Vector jx_carom(zero_vector.GetData(), zero_vector.Size(), true, false);  
    // CAROM::Vector jy_carom(zero_vector.GetData(), zero_vector.Size(), true, false);

    std::unique_ptr<CAROM::Vector> reducedRHS = spatialbasis->transposeMult(B_carom); // Reduced version of B : phi^T*B

    CAROM::Vector reducedHn(numColumnRB, false);  // Reduced version of Hn
    reducedHn = 0;
    CAROM::Vector reducedHnm1(numColumnRB, false);
    reducedHnm1 = 0;

    CAROM::Vector reducedBn(numColumnRB, false);  // Reduced version of Bn 
    reducedBn = 0;
    CAROM::Vector reducedBnm1(numColumnRB, false);
    reducedBnm1 = 0;

    CAROM::Vector reducedBdtn(numColumnRB, false);
    reducedBdtn = 0;

    // CAROM::Vector reducedEx(numColumnRB, false);
    // reducedEx = 0;
    // CAROM::Vector reducedEy(numColumnRB, false);
    // reducedEy = 0;

    // CAROM::Vector reducedExnm1(numColumnRB, false);
    // reducedExnm1 = 0;
    // CAROM::Vector reducedEynm1(numColumnRB, false);
    // reducedEynm1 = 0;

    // CAROM::Vector reducedJxnm1(numColumnRB, false);
    // reducedJxnm1 = 0;
    // CAROM::Vector reducedJynm1(numColumnRB, false);
    // reducedJynm1 = 0;




    assembleTimer.Stop();

    std::cout << "Number of dof : " << fespace->GetNDofs() << std::endl;
    std::cout << "VSize : " << fespace->GetTrueVSize() << std::endl;

 
    socketstream sout;

    char vishost[] = "localhost";
    int visport = 19916;

    if (visualization)  // Initialize the ploting of the results
    {
        sout.open(vishost, visport);
        sout << "parallel " << num_procs << " " << myid << "\n";
        sout.precision(8);
        sout << "solution\n" << pmesh << Jn << "\nkeys j\n" << std::flush;

    }


    int vis_steps = 10;  // Number of iterations between two visualizations
    for (int step = 0; step < num_steps; step++)
    {   
        iter++;

        t += Ts;                    // New time 
        bdr_coeff.SetTime(t);       // Update boundary conditions fonction's time parameters
        Hn.ProjectBdrCoefficient(bdr_coeff, dir_bdr); // Setting new boundary conditions on Hn

        b = 0;
        b.Assemble();

        r1.FormLinearSystem(ess_tdof_list, Hn, b, A, X, B);  // On doit Recalculer rhs à la fin de la boucle
        // Here B contains the bondary conditions only, it'll be used to add thoses boundary conditions to rhs later after being reduced
        CAROM::Vector B_carom_bdr(B.GetData(), B.Size(), true, false); 
        std::unique_ptr<CAROM::Vector> B_bdr_reduced = spatialbasis->transposeMult(B_carom_bdr);
    
        // reducedRHS = reducedR2 * reduced_hnm1 + reducedM * reducedDBdT
        std::unique_ptr<CAROM::Vector> reduced_R2Hnm1 = ReducedR2.mult(reducedHnm1);
        std::unique_ptr<CAROM::Vector> reduced_MdBdt = ReducedM.mult(reducedBdtn);
        *reducedRHS = *reduced_R2Hnm1;
        *reducedRHS += *reduced_MdBdt;

        *reducedRHS += *B_bdr_reduced;

        solveTimer.Start();
        invReducedA.mult(*reducedRHS, reducedHn);

        solveTimer.Stop();

        // Update of db_dt with : (db_dt)_n = B1*Hn + B2*Hn-1 + B3*(db_dt)_n-1
        reducedBdtn *= B3;
        reducedHnm1 *= B2;
        reducedHn *= B1;

        reducedBdtn += reducedHnm1;
        reducedBdtn += reducedHn;

        reducedHn *= 1/B1;
        reducedHnm1 *= 1/B2;

        // Bn = C1 Hn + C2 Hnm1 + C3 Bnm1
        reducedHn *= C1;
        reducedHnm1 *= C2;
        reducedBnm1 *= C3;

        reducedBn = reducedHn;
        reducedBn += reducedHnm1;
        reducedBn += reducedBnm1; // reduecedBn up to date

        reducedHn *= 1/C1;
        reducedHnm1 = reducedHn;
        reducedBnm1 = reducedBn;


        // std::unique_ptr<CAROM::Vector> reduced_Jx = ReducedCx.mult(reducedHn);
        // std::unique_ptr<CAROM::Vector> reduced_Jy = ReducedCy.mult(reducedHn);

        // reducedExnm1 *= A3;
        // reducedJxnm1 *= A2;
        // *reduced_Jx *= A1;
        // reducedEx = *reduced_Jx;
        // reducedEx += reducedJxnm1;
        // reducedEx += reducedExnm1;

        // reducedEynm1 *= A3;
        // reducedJynm1 *= A2;
        // *reduced_Jy *= A1;
        // reducedEy = *reduced_Jy;
        // reducedEy += reducedJynm1;
        // reducedEy += reducedEynm1;

        // *reduced_Jx *= 1/A1;
        // *reduced_Jy *= 1/A1;
        
        // reducedJxnm1 = *reduced_Jx;
        // reducedJynm1 = *reduced_Jy;
        // reducedExnm1 = reducedEx;
        // reducedEynm1 = reducedEy;

        /*
        Jn_temp = Jn;
        Jn_temp *=A1;
        Jnm1 *= A2;
        Enm1 *= A3;
        En = 0;
        En+=Jn_temp;
        En+=Jnm1;
        En+=Enm1;
        */

        spatialbasis->mult(reducedHn, X_carom);  // Getting full version of Hn in X_carom (compatible libRom)
        spatialbasis->mult(reducedBn, Bn_carom);  // Getting full version of Hn in X_carom (compatible libRom)
        // spatialbasisEx->mult(*reduced_Jx, jx_carom);
        // spatialbasisEy->mult(*reduced_Jy, jy_carom);
        // spatialbasisEx->mult(reducedEx, Ex_carom);
        // spatialbasisEy->mult(reducedEy, Ey_carom);

        for (int i = 0; i<Hn.Size(); i++){       
            Hn(i) = X_carom(i);                  // Getting Hn in a mfem::Vector
            Bn(i) = Bn_carom(i);
            // Pn(i) = jx_carom(i) * Ex_carom(i) + jy_carom(i) * Ey_carom(i);
            // Jnx(i) = jx_carom(i);
        }


        Hn.ProjectBdrCoefficient(bdr_coeff, dir_bdr);  // Ensuring bdr conditions (may not be necessary)

        // // ***** computation of magnetic flux ***********
        // GridFunctionCoefficient B_coeff(&Bn);
        // GridFunction one(fespace);
        // one = 1;
        // LinearForm lf_flux(fespace);
        // lf_flux.AddDomainIntegrator(new DomainLFIntegrator(B_coeff));
        // lf_flux.Assemble();

        // real_t flux = lf_flux(one);

        // GridFunctionCoefficient P_eddy_coeff(&Pn);
        // LinearForm lf(fespace);
        // lf.AddDomainIntegrator(new DomainLFIntegrator(P_eddy_coeff));
        // lf.Assemble();
        // real_t P_eddy_tot = 2*M_PI*lf(one);
        // real_t vol = M_PI * height * (Rout*Rout - Ri*Ri);
        // real_t P_eddy = P_eddy_tot/vol;

        real_t flux = 0;
        real_t P_eddy = 0;
        real_t flux_test = 0;
        

        // Cx.Mult(Hn, Jnx);
        // Cy.Mult(Hn, Jny);
        // for (int i = 0; i<Hn.Size(); i++) {
        //     Jn(i) = Jnx(i);
        //     Jn(i + Hn.Size()) = Jny(i);
        // }

        curl.Mult(Hn, Jn);

        data_file << t << ";" << P_eddy << ";" << flux << ";" << NI_func(t) <<";" << flux_test << std::endl;




        if (step % vis_steps == 0 )
        {       
            std::cout << "Time step: " << step << ", Time: " << t << std::endl;
            string name = "data/reduced/Hn_reduced_"+to_string(step);  
            X_carom.print(name.c_str());
            // Vector lf_vect = lf;
            // for (int i = 0; i<50; i++){
            //     std::cout << lf_carom(i) << std::endl;
            // }
            if (visualization)
            {
                // Cy.Mult(Hn, Jnx_test);
                // std::cout << Cx.Height() << std::endl;
                sout.precision(8);
                sout << "solution\n"
                     << pmesh << Jn
                     << "window_title 'Champ H'"
                     << std::flush;
                std::cin.get();
            }
            
        }
    }
}



PowerLossCoefficient_TD::PowerLossCoefficient_TD(const FiniteElementSpace *fespace_, CurlCustomCoefficient &J_, VectorGridFunctionCoefficient &E_)
    : fespace(fespace_),
      J(J_), J_vect(fespace_->GetMesh()->SpaceDimension()),
      E(E_), E_vect(fespace_->GetMesh()->SpaceDimension()) {}

real_t PowerLossCoefficient_TD::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    Vector x;               // Vector coordinates
    T.Transform(ip, x);     // Get the global coordinates in vector x from integration point's coordinates in the element referential
    J.Eval(J_vect, T, ip);    // Get from J_r (Coefficient) the value at the point ip in J_r_vect
    E.Eval(E_vect, T, ip);    // same for E
    return  (J_vect*E_vect) * x[0];  // J.E * r (Cylindrical coordinates)
}




// void ComputeCurl(const mfem::HypreParMatrix& Mass,
//                            const mfem::HypreParMatrix& Kx,
//                            CAROM::Matrix& C) // Output: dense, local
// {
//     const int n = Mass.NumRows();
//     const int m = Kx.NumCols();

//     MFEM_VERIFY(Kx.NumRows() == n && C.numRows() == n && C.numColumns() == m,
//                 "Dimensions incompatibles entre Mass, Kx et C.");

//     // Setup solveur
//     mfem::CGSolver solver(MPI_COMM_WORLD);
//     mfem::HypreBoomerAMG amg(Mass); // Préconditionneur
//     solver.SetOperator(Mass);
//     solver.SetPreconditioner(amg);
//     solver.SetRelTol(1e-10);
//     solver.SetMaxIter(500);
//     solver.SetPrintLevel(0);

//     // Vecteurs de travail
//     mfem::Vector ej(m);      // base canonique
//     mfem::Vector Kxej(n);    // Kx * ej
//     mfem::Vector xj(n);      // solution M^{-1} * Kxej

//     for (int j = 0; j < m; ++j)
//     {
//         ej = 0.0;
//         ej[j] = 1.0;

//         Kx.Mult(ej, Kxej);     // j-ième colonne de Kx
//         solver.Mult(Kxej, xj); // xj = M^{-1} * Kxej

//         // Stockage dans C
//         for (int i = 0; i < n; ++i)
//         {
//             C(i, j) = xj[i];
//         }
//     }
// }
