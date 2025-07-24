#include "solver_TD_FOM.hpp"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"
#include <iostream>
#include <cmath>
#include <chrono>


using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

void CustomCurlInterpolatorX::AssembleElementMatrix2(const FiniteElement &h1_fe,
                            const FiniteElement &out_fe,
                            ElementTransformation &Trans,
                            DenseMatrix &elmat)
{
    const int dof = h1_fe.GetDof();         // Number of degrees of freedom inside the element with the input formulation (Here H1)
    const int dof_out = out_fe.GetDof();    // Number of degrees of freedom inside the element with the output formulation (Here H1 too)
    const int dim = out_fe.GetDim();        // Here dim is 2
    DenseMatrix dshape(dof, dim);           // Matrix that will contain the gradient of each dofs. Size : dofs x dim
    DenseMatrix dshapedxt(dof, dim);        // Same but in the physical domain. Size : dofs x dim
    DenseMatrix invdfdx(dim, dim);          // Adjugate of the Jacobian of the change of variable between physical and element space
    
    elmat.SetSize(dof_out, dof);            // Local matrix of the contribution of the degrees of freedom. Size : dof_out x dof

    for (int i = 0; i < dof_out; i++)       // Looping the gradient of the shape functions at the coordinates of the dof i  
    {
        const IntegrationPoint &ip = out_fe.GetNodes().IntPoint(i);
        Trans.SetIntPoint(&ip);

        Vector phys(Trans.GetSpaceDim());      // (x,y) in the element's reference
        Trans.Transform(ip, phys);             // phys = (r,z) of the physical space 
        h1_fe.CalcDShape(ip, dshape);          // Gradient of all shape function of the element in the element's space
        CalcAdjugate(Trans.Jacobian(), invdfdx); 
        Mult(dshape, invdfdx, dshapedxt);      // Computing the gradient in the physical space

        real_t w = 1.;      
        // H1 space (not weighted by the volume of the element), thus it is necessary to remove the weihting of the volume of the element
        if (out_fe.GetMapType() == FiniteElement::VALUE)  
        {
            Trans.SetIntPoint(&ip); 
            w /= Trans.Weight();  
        }

        for (int j = 0; j < dof; j++)   
        {
            // If the dof has not been computed yet by an other element, we apply the value -dh/dz, else 0, to avoid cumulative effect
            // A better approach would be to compute the average of the contributions of the differents elements on one dof and not only one
            elmat(i, j) = (set.find(phys(0)+ 10000* phys(1)) != set.end() ? 0 : -dshapedxt(j, 1) * w);
        }
        // Very simple and naive "hash function". Could be improved to guarantee the unicity
        // Here it is working as the sizes of the elements are "big enough" 
        set.insert(phys(0)+ 10000* phys(1));  // Add the 'dof' to the set
    }
}


void CustomCurlInterpolatorY::AssembleElementMatrix2(const FiniteElement &h1_fe,
                            const FiniteElement &out_fe,
                            ElementTransformation &Trans,
                            DenseMatrix &elmat)
{
    const int dof = h1_fe.GetDof();         // Number of degrees of freedom inside the element with the input formulation (Here H1)
    const int dof_out = out_fe.GetDof();    // Number of degrees of freedom inside the element with the output formulation (Here H1 too)
    const int dim = out_fe.GetDim();        // Here dim is 2
    DenseMatrix dshape(dof, dim);           // Matrix that will contain the gradient of each dofs. Size : dofs x dim
    DenseMatrix dshapedxt(dof, dim);        // Same but in the physical domain. Size : dofs x dim
    DenseMatrix invdfdx(dim, dim);          // Adjugate of the Jacobian of the change of variable between physical and element space
    
    elmat.SetSize(dof_out, dof);            // Local matrix of the contribution of the degrees of freedom

    for (int i = 0; i < dof_out; i++)       // Looping on every degree of freedom of the element
    {
        const IntegrationPoint &ip = h1_fe.GetNodes().IntPoint(i);
        Trans.SetIntPoint(&ip);

        Vector phys(Trans.GetSpaceDim());      // (x,y) in the element's reference
        Trans.Transform(ip, phys);             // phys = (r,z) of the physical space 
        h1_fe.CalcDShape(ip, dshape);          // Gradient of all shape function of the element in the element's space
        CalcAdjugate(Trans.Jacobian(), invdfdx); 
        Mult(dshape, invdfdx, dshapedxt);      // Computing the gradient in the physical space

        real_t w = 1.; 
        // H1 space (not weighted by the volume of the element), thus it is necessary to remove the weihting of the volume of the element
        if (out_fe.GetMapType() == FiniteElement::VALUE)
        {
            Trans.SetIntPoint(&ip);
            w /= Trans.Weight();
        }
        for (int j = 0; j < dof; j++)
        {
            // If the dof has not been computed yet by an other element, we apply the value -dh/dz, else 0, to avoid cumulative effect
            // A better approach would be to compute the average of the contributions of the differents elements on one dof and not only one
            elmat(i, j) = (set.find(phys(0)+ 10000* phys(1)) != set.end() ? 0 : (i==j ? dshapedxt(j, 0) * w + 1/phys(0) : dshapedxt(j, 0) * w));
        }
        // Very simple and naive "hash function". Could be improved to guarantee the unicity
        // Here it is working as the sizes of the elements are "big enough"
        set.insert(phys(0)+ 10000* phys(1)); // Add the 'dof' to the set
    }
}

// Inspired from ComputeCtAB of librom : https://librom.readthedocs.io/en/latest/
void ComputeCtAB_Sparse(const mfem::SparseMatrix& A,
                        const CAROM::Matrix& B,  // Distributed matrix
                        const CAROM::Matrix& C,  // Distributed matrix
                        CAROM::Matrix& CtAB)     // Non-distributed (local) matrix
{
    MFEM_VERIFY(B.distributed() && C.distributed() && !CtAB.distributed(),
                "In ComputeCtAB, B and C must be distributed, but not CtAB.");

    const int num_rows = B.numRows();       // Local rows of B
    const int num_cols = B.numColumns();    // Number of columns in B
    const int num_rows_A = A.Height();      // Number of rows in A

    MFEM_VERIFY(C.numRows() == num_rows_A, "C rows must match A rows.");

    mfem::Vector Bvec(num_rows);
    mfem::Vector ABvec(num_rows_A);

    CAROM::Matrix AB(num_rows_A, num_cols, true);

    for (int i = 0; i < num_cols; ++i) {
        for (int j = 0; j < num_rows; ++j) {
            Bvec[j] = B(j, i);
        }
        A.Mult(Bvec, ABvec);
        for (int j = 0; j < num_rows_A; ++j) {
            AB(j, i) = ABvec[j];
        }
    }
    C.transposeMult(AB, CtAB);
}



void TD_sim_offline(ParMesh &pmesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization, bool save_snapshot, int num_procs, int myid){

    // ****************** Parameters definitions *********************
    real_t Ts = t_f/(num_steps);    // Time step

    // ******* geometric parameters ******          Representation of the Mesh
                                    //                      ____________
    real_t Ri = 9.6e-3/2.0;         //           ^         |            |
    real_t height = 7.59e-3;        //           |         |            |            
    real_t w = 5.3e-3;              //           |         |            |            
    real_t Rout = Ri + w;           //           | Height  |            |
    real_t Rm = (Rout - Ri) / 2;    //           |         |            |
    int iter = 0;                   //           |         |____________|
                                    //                           w
                                    //      Ri -----------><------------>
                                    //    Rout ------------------------->

    // ****** Ferrite parameters ******
    //  Measured parameters :
    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;
    real_t mu = ferrite.mu; 
    
    // Fitted parameter : 
    real_t tau = 1/(2 * M_PI * 1.8e6);
    // *******************************


    // *********** Time Domain pararmeters *********
    real_t t = 0.0;

    // Z transform of rho_eq -> (A1 + A2.z^(-1)) / (1 - A3.z^(-1))
    real_t A1 = (rho * (sigma * Ts + 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A2 = (rho * (sigma * Ts - 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A3 = -(Ts*sigma -2*eps) / (2*eps + Ts*sigma);

    // Z transform of s.mu_eq -> (B1 + B2.z^(-1)) / (1 - B3.z^(-1)) with s the Laplace variable
    real_t B1 = 2*mu / (Ts + 2*tau);
    real_t B2 = -(2*mu) / (Ts + 2*tau);
    real_t B3 = -(Ts - 2*tau) / (Ts + 2*tau);

    // Z transform of  mu_eq : (C1 + C2.z^(-1)) / (1 - C3.z^(-1))
    real_t C1 = Ts*mu/(Ts+2*tau);
    real_t C2 = Ts*mu/(Ts+2*tau);
    real_t C3 = -(Ts-2*tau)/(Ts+2*tau);


    int order = 2;                      // elements order : 2
    int dim = pmesh.Dimension();         // Mesh dimension : Here 2
    int id = 0;                         // id of the main mpi rank

    // Boundary conditions function
    auto bdr_func = [&](const Vector &x, real_t t) -> real_t
    {
        real_t r = x(0); // Accès au premier élément du vecteur x
        return NI_func(t) / (2 * M_PI * r);
    };

    // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];
    };

    // Function for the correcting factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };
    
    // Using a H1 finite element collection
    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(&pmesh, fec);        // Finite elements space for the magnetic field (Theta component)
    ParFiniteElementSpace *fespace_E = new ParFiniteElementSpace(&pmesh, fec, dim); // Finite elements space fot the electrical field (r and z components)
    
    HYPRE_Int size = fespace->GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of finite element unknowns: " << size << endl;
    }

    // Initialization of the differents gridFunctions
    ParGridFunction Hn(fespace);        // Magnetic field at time step n
    ParGridFunction Hnm1(fespace);      // Magnetic field at time step n-1
    ParGridFunction db_dt(fespace);     // dB/dt at time step n
    ParGridFunction db_dtm1(fespace);   // dB/dt at time step n-1

    ParGridFunction Bn(fespace);        // Magnetic flux density field at time step n
    ParGridFunction Bnm1(fespace);      // Magnetic flux density field at time step n-1

    ParGridFunction En(fespace_E);      // Electrical field at time step n
    ParGridFunction Enm1(fespace_E);    // Electrical field at time step n-1

    ParGridFunction Jn(fespace_E);      // Current density at time step n
    ParGridFunction Jn_x(fespace);      // r component of the current density (used as snapshots for the reduced model)
    ParGridFunction Jn_y(fespace);      // z component of the current density (used as snapshots for the reduced model)
    ParGridFunction Jnm1(fespace_E);    // Current density at time step n-1
    
    ParGridFunction Pn(fespace);        // Power density (En.Jn)

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
    dir_bdr = 1; // All borders have Dirichlet boundary conditions (dbc)
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // Getting the list of degrees of freedom with dbc

    FunctionCoefficient bdr_coeff([&](const Vector &x) { return bdr_func(x, t); });
    // Hnm1.ProjectBdrCoefficient(bdr_coeff, dir_bdr); // Project the bdc on the previous time step (not necessary as the source term is 0 before) 

    // **** Coefficients for bilinear forms ***********
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

    // Weak formulation : r1(h,h') = r2(hnm1,h') + m(dBn-1/dt, h') 

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
    

    // ******** Parameters for the linear system *********
    Vector rhs(fespace->GetTrueVSize());
    Vector R2Hnm1(fespace->GetTrueVSize());
    Vector MdBdt(fespace->GetTrueVSize());

    Vector x(fespace->GetTrueVSize()); // For defining the linear system, considering Dirichlet's bdc

    // Using Conjugate Gradient Solver
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


    // ********* libROM parameters (For saving the snapshots for the reduced model) ****************
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    CAROM::Options* options;
    int numRowRB, numColumnRB;

    // Basis for the reduced magnetic field 
    const std::string basisName = "data/basis";
    const std::string basisFileName = basisName + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasis;
    CAROM::BasisGenerator *generator;
    
    // Basis for the reduced electrical field 
    const std::string basisNameEx = "data/basisEx";
    const std::string basisNameEy = "data/basisEy";
    const std::string basisFileNameEx = basisNameEx + std::to_string(id);
    const std::string basisFileNameEy = basisNameEy + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasisEx;
    std::unique_ptr<const CAROM::Matrix> spatialbasisEy;
    CAROM::BasisGenerator *generatorEx;
    CAROM::BasisGenerator *generatorEy;

    options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots, update_right_SV);

    // generator objects that will generate the basis from the snapshots using Singular Value Decomposition (SVD)   
    generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
    generatorEx = new CAROM::BasisGenerator(*options, isIncremental, basisFileNameEx);
    generatorEy = new CAROM::BasisGenerator(*options, isIncremental, basisFileNameEy);


    // Linear Form the magnetic flux : lf_flux^T . Bn = intergral of Bn = phi_n
    ConstantCoefficient ones(1.);
    LinearForm lf_flux(fespace);
    lf_flux.AddDomainIntegrator(new DomainLFIntegrator(ones));
    lf_flux.Assemble();

    // Linear Form the power losses: lf_p^T . Pn = integral of Pn
    LinearForm lf_p(fespace);
    lf_p.AddDomainIntegrator(new DomainLFIntegrator(ones));
    lf_p.Assemble();

    std::cout << "Starting FOM Simulation : " << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();  // Time of the begining of the simulation
    int vis_steps = 100;                                   // Number of iterations between two visualizations (if visualization)
    // ***************  Time iterations *****************
    for (int step = 0; step < num_steps; step++)
    {   
        iter++;
        t += Ts;                    // Update time
        bdr_coeff.SetTime(t);       // Update boundary conditions fonction's time parameters
        Hn.ProjectBdrCoefficient(bdr_coeff, dir_bdr); // Setting new boundary conditions on H

        // Mult M * dBdt → MdBdt
        m.Mult(db_dtm1, MdBdt);

        // Mult R2 * Hnm1 → R2Hnm1
        r2.Mult(Hnm1, R2Hnm1);

        // RHS = R2Hnm1 + MdBdt = R2 . Hnm1 + M . dBnm1/dt
        rhs = R2Hnm1;
        rhs += MdBdt;
        
        HypreParMatrix A_sys;
        Vector X, B;

        // Form linear system taking boundary conditions into account
        // Here A_sys is a modified version of r1 taking the dbc into account
        // Here A_sys is a modified version of rhs taking the dbc into account
        r1.FormLinearSystem(ess_tdof_list, Hn, rhs, A_sys, X, B);
        

        // Solving Linear system : A_sys * X = B
        solver.SetOperator(A_sys);
        solver.Mult(B, X);
        r1.RecoverFEMSolution(X, rhs, Hn);  // Export the solution of the system into the GridFunction Hn


        // *** Update of the other GridFunctions from the solution Hn ***

        // Update of db_dt with : (db_dt)_n = B1*Hn + B2*Hn-1 + B3*(db_dt)_n-1  <-> db/dt = s.mu_eq Hn (see Z transform of s.mu_eq above)
        db_dtm1 *= B3;
        ParGridFunction H_temp(fespace); // temporary GridFunction to avoid modifying Hn
        H_temp = Hn;
        H_temp *= B1;  // Hn*B1

        Hnm1 *= B2;     // Hn-1 * B2

        db_dt = H_temp;
        db_dt += Hnm1;
        db_dt += db_dtm1;  // db_dt up to date

        Hnm1 /= B2;         // Restoring Hnm1 and db_dtm1 for the update of Bn
        db_dtm1 /= B3;
        H_temp /= B1;
    
        // ***** Update of Bn ******  
        // Bn = C1 Hn + C2 Hnm1 + C3 Bnm1   <->  B = mu_eq Hn (see Z transform of mu_eq)
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
        real_t flux = lf_flux(Bn);

        // *********** Computation of J and E *********** 
        CurlCustomCoefficient J_coeff(&Hn);
        Jn.ProjectCoefficient(J_coeff);     // Jn = curl(Hn)

        // En = A1 . Jn + A2 . Jn-1 + A3 . En-1
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

        VectorGridFunctionCoefficient E_coeff(&En);  // Coefficient used for the power losses
        PowerLossCoefficient_TD P_eddy_coeff(fespace_E, J_coeff, E_coeff);  // Coefficient of the power losses density
        Pn.ProjectCoefficient(P_eddy_coeff);  // Compute the power losses density in the gridFunction Pn

        real_t P_eddy_tot = 2*M_PI*lf_p(Pn);  // Losses on the volume
        real_t vol = M_PI * height * (Rout*Rout - Ri*Ri);
        real_t P_eddy = P_eddy_tot/vol;       // Average power density

        // Writing the results of the power losses and the flux
        data_file << t << ";" << P_eddy << ";" << flux << ";" << NI_func(t) << std::endl;

        // Actions each vis_step (visualization or saving snapshots)
        if (step % vis_steps == 0 )
        {   
            std::cout << "Time step: " << step << ", Time: " << t << std::endl;
            if (save_snapshot)
            {
                // For the reduced model, we used the each component of the current density as snaphots
                for (int i = 0; i<Hn.Size(); i++) {
                    Jn_x(i) = Jn(i);
                    Jn_y(i) = Jn(i+Hn.Size());
                }
                
                bool addSample = generator->takeSample(Hn.GetData());       // Add the solution Hn in the snapshots
                bool addSampleEx = generatorEx->takeSample(Jn_x.GetData()); // Add the solution Jn in the snapshots (each components)
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
                sout_b << "solution\n" << pmesh << Pn <<  "window_title 'Field B'" << std::flush;

                // std::cin.get(); // Uncomment for pausing the visualization at each step
            }
        }
        
    }
    auto t1 = std::chrono::high_resolution_clock::now();  // Capture the time after the simulation
    std::chrono::duration<double> dt = t1 - t0;           // Get the duration of the simulation
    std::cout << "Elapsed time for solving FOM: " <<  dt.count() << " s\n";
    std::cout << "FOM Simulation has ended successfully." << std::endl;
    // Save the snaphots and generate the basis using SVD
    if (save_snapshot) {
        std::cout << "Generating snapshots files : " << std::endl;
        generator->writeSnapshot();
        generatorEx->writeSnapshot();
        generatorEy->writeSnapshot();
        generator->endSamples();
        generatorEx->endSamples();
        generatorEy->endSamples();
    }

    data_file.close();
    // Free memory
    delete generator;
    delete generatorEx;
    delete generatorEy;
    delete options;
    delete fespace_E;
    delete fespace;
    delete fec;
}

void TD_sim_online(ParMesh &pmesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization, int num_procs, int myid){
    
    // ****************** Parameters definitions *********************
    real_t Ts = t_f/(num_steps);    // Time step

    // ******* geometric parameters ******          Representation of the Mesh
                                    //                      ____________
    real_t Ri = 9.6e-3/2.0;         //           ^         |            |
    real_t height = 7.59e-3;        //           |         |            |            
    real_t w = 5.3e-3;              //           |         |            |            
    real_t Rout = Ri + w;           //           | Height  |            |
    real_t Rm = (Rout - Ri) / 2;    //           |         |            |
    int iter = 0;                   //           |         |____________|
                                    //                           w
                                    //      Ri -----------><------------>
                                    //    Rout ------------------------->

    // ****** Ferrite parameters ******
    //  Measured parameters :
    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;
    real_t mu = ferrite.mu; 
    
    // Fitted parameter : 
    real_t tau = 1/(2 * M_PI * 1.8e6);
    // *******************************


    // *********** Time Domain pararmeters *********
    real_t t = 0.0;

    // Z transform of rho_eq -> (A1 + A2.z^(-1)) / (1 - A3.z^(-1))
    real_t A1 = (rho * (sigma * Ts + 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A2 = (rho * (sigma * Ts - 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A3 = -(Ts*sigma -2*eps) / (2*eps + Ts*sigma);

    // Z transform of s.mu_eq -> (B1 + B2.z^(-1)) / (1 - B3.z^(-1)) with s the Laplace variable
    real_t B1 = 2*mu / (Ts + 2*tau);
    real_t B2 = -(2*mu) / (Ts + 2*tau);
    real_t B3 = -(Ts - 2*tau) / (Ts + 2*tau);

    // Z transform of  mu_eq : (C1 + C2.z^(-1)) / (1 - C3.z^(-1))
    real_t C1 = Ts*mu/(Ts+2*tau);
    real_t C2 = Ts*mu/(Ts+2*tau);
    real_t C3 = -(Ts-2*tau)/(Ts+2*tau);


    int order = 2;                      // elements order : 2
    int dim = pmesh.Dimension();         // Mesh dimension : Here 2
    int id = 0;                         // id of the main mpi rank


    // Boundary conditions function
    auto inv_r_func = [&](const Vector &x, real_t t) -> real_t
    {
        return 1. / (2 * M_PI * x[0]);
    };

    // // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];
    };

    // Function for the correcting factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };
    


   // Using a H1 finite element collection
    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(&pmesh, fec);        // Finite elements space for the magnetic field (Theta component)
    ParFiniteElementSpace *fespace_E = new ParFiniteElementSpace(&pmesh, fec, dim); // Finite elements space fot the electrical field (r and z components)
    
    HYPRE_Int size = fespace->GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of finite element unknowns: " << size << endl;
    }

    // Initialization of the differents gridFunctions
    ParGridFunction Hn(fespace);        // Magnetic field at time step n
    ParGridFunction Hnm1(fespace);      // Magnetic field at time step n-1
    ParGridFunction db_dt(fespace);     // dB/dt at time step n
    ParGridFunction db_dtm1(fespace);   // dB/dt at time step n-1

    ParGridFunction Bn(fespace);        // Magnetic flux density field at time step n
    ParGridFunction Bnm1(fespace);      // Magnetic flux density field at time step n-1

    ParGridFunction En(fespace_E);      // Electrical field at time step n
    ParGridFunction Enx(fespace);       // r component of the electrical field (used as snapshots for the reduced model)
    ParGridFunction Eny(fespace);       // z component of the electrical field (used as snapshots for the reduced model)
    ParGridFunction Enm1(fespace_E);    // Electrical field at time step n-1

    ParGridFunction Jn(fespace_E);      // Current density at time step n
    ParGridFunction Jnx(fespace);       // r component of the current density (used as snapshots for the reduced model)
    ParGridFunction Jny(fespace);       // z component of the current density (used as snapshots for the reduced model)
    ParGridFunction Jnm1(fespace_E);    // Current density at time step n-1
    
    ParGridFunction Pn(fespace);        // Power density (En.Jn)

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

    // Definition of the Boundary conditions :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(pmesh.bdr_attributes.Max());
    dir_bdr = 1; // All borders have Dirichlet boundary conditions (dbc)
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // Getting the list of degrees of freedom with dbc

    // FunctionCoefficient bdr_coeff([&](const Vector &x) { return bdr_func(x, t); });
    ConstantCoefficient zero_coeff(0.);
    // Hnm1.ProjectBdrCoefficient(zero_coeff, dir_bdr);  // Not necessary as the initial bdc is 0

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

    // Weak formulation : r1(h,h') = r2(hnm1,h') + m(dBn-1/dt, h') 

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

    // Initialization of Hypre matrices 
    HypreParMatrix A, R2, M;

    ParLinearForm b(fespace);
    b = 0;
    b.Assemble();
    Vector X, B;
    Vector zero_vector(Hn.Size());
    
    Vector rhs(fespace->GetTrueVSize());
    Vector R2Hnm1(fespace->GetTrueVSize());
    Vector MdBdt(fespace->GetTrueVSize());

    // Applying 1/(2*pi*r) on the border. Will be multiplied by NI later (Directly in the reduced basis as it's linear)
    FunctionCoefficient inv_r_coeff(inv_r_func);
    Hn.ProjectBdrCoefficient(inv_r_coeff, dir_bdr);


    Array<int> ess_tdof_list_fake; // Used to generate the Matrixes of the bilinear forms r2 and m 

    r1.FormLinearSystem(ess_tdof_list, Hn, b, A, X, B); // Here using FormLinearSystem to initialize A and B, such that we have AX = B
                                                        // B contains the dbc applied, which are 1/(2*pi*r)
                                                        // but the real dbc are : NI(t)/(2*pi*r). However, as the 
                                                        // problem is linear, we can project the bdc (here B) into the reduced space
                                                        // and then multiply it by NI. Which costs less ressources
    
    r2.FormSystemMatrix(ess_tdof_list_fake, R2); // Generating the matrixes R2 and M from the bilinear forms
    m.FormSystemMatrix(ess_tdof_list_fake, M);   // --------------------------------


    // ---------- Generation of curl matrices -----------------
    // Generation of Cx matrix using a customed interpolator.
    // The matrix Cx is defined such that Cx.Hn = Jnx. 
    // It is necessary to use such a matrix for being able to project 
    // this linear operation in the reduced basis
    ParDiscreteLinearOperator curlx(fespace, fespace);
    curlx.AddDomainInterpolator(new CustomCurlInterpolatorX());
    curlx.Assemble();
    curlx.Finalize();
    HypreParMatrix *Cx = curlx.ParallelAssemble();
    
    ParDiscreteLinearOperator curly(fespace, fespace);
    curly.AddDomainInterpolator(new CustomCurlInterpolatorY());
    curly.Assemble();
    curly.Finalize();
    HypreParMatrix *Cy = curly.ParallelAssemble();
    
    // ----------------------------------------------------------

    
    // File for saving the values of the power and flux for each iterations
    std::string name = "./data/reduced/TD_" + std::to_string(1) + ".csv";   
    std::ofstream data_file(name);                         
    data_file << "t;p_eddy;flux;NI\n0;0;0;0\n";  // Intialising the file with coluns names and first values to 0


    // ********* libROM parameters (For saving the snapshots for the reduced model) ****************
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    CAROM::Options* options;
    int numRowRB, numColumnRB;

    // Basis for the reduced magnetic field 
    const std::string basisName = "data/basis";
    const std::string basisFileName = basisName + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasis;
    
    // Basis for the reduced electrical field 
    const std::string basisNameEx = "data/basisEx";
    const std::string basisNameEy = "data/basisEy";
    const std::string basisFileNameEx = basisNameEx + std::to_string(id);
    const std::string basisFileNameEy = basisNameEy + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasisEx;
    std::unique_ptr<const CAROM::Matrix> spatialbasisEy;

    // Reading the bases for projecting in the reduced space 
    CAROM::BasisReader reader(basisFileName);      // Get Snapshot in a BasisReader object
    CAROM::BasisReader readerEx(basisFileNameEx);  // Get Snapshot in a BasisReader object
    CAROM::BasisReader readerEy(basisFileNameEy);  // Get Snapshot in a BasisReader object

    spatialbasis = reader.getSpatialBasis();     // Get the reduced basis of the magnetic field from the snapshots
    spatialbasisEx = readerEx.getSpatialBasis(); // Get the reduced basis of Ex from the snapshots
    spatialbasisEy = readerEy.getSpatialBasis(); // Get the reduced basis of Ey from the snapshots

    numRowRB = spatialbasis->numRows();          // Number of dof of the classic FEM model
    numColumnRB = spatialbasis->numColumns();    // Number of reduced dof, here equal to the number of snapshots
    if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB, numColumnRB);

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
    CAROM::Matrix ReducedCx(numColumnRB, numColumnRB, false);
    CAROM::Matrix ReducedCy(numColumnRB, numColumnRB, false);


    ComputeCtAB(A, *spatialbasis, *spatialbasis, ReducedA);    // Reduced A = phi^t * A * phi
    ComputeCtAB(A, *spatialbasis, *spatialbasis, invReducedA); // inverse of Reduced A 
    invReducedA.inverse();
    ComputeCtAB(R2, *spatialbasis, *spatialbasis, ReducedR2); // Reduced R2
    ComputeCtAB(M, *spatialbasis, *spatialbasis, ReducedM);   // Reduced M
    invReducedA.print("data/A_inv");
    ReducedA.print("data/A");

    ComputeCtAB(*Cx, *spatialbasis, *spatialbasisEx, ReducedCx);    // Reduced Cx
    ComputeCtAB(*Cy, *spatialbasis, *spatialbasisEy, ReducedCy);    // Reduced Cy

    CAROM::Vector B_carom(B.GetData(), B.Size(), true, false);     // B Compatible with libRom, share data with mfem vector
    CAROM::Vector X_carom(X.GetData(), X.Size(), true, false);     // X Compatible with libRom 
    CAROM::Vector Bn_carom(Bn.GetData(), Bn.Size(), true, false);  // X Compatible with libRom
    
    CAROM::Vector Ex_carom(Enx.GetData(), Enx.Size(), true, false);  // Ex compatible wit libRom
    CAROM::Vector Ey_carom(Eny.GetData(), Eny.Size(), true, false);  // Ey compatible wit libRom

    CAROM::Vector jx_carom(Jnx.GetData(), Jnx.Size(), true, false);  // Jx compatible wit libRom 
    CAROM::Vector jy_carom(Jny.GetData(), Jny.Size(), true, false);  // Jy compatible wit libRom


    // Initialization of each "GridFunction/Vector" in the reduced space 
    std::unique_ptr<CAROM::Vector> reducedRHS = spatialbasis->transposeMult(B_carom); // Reduced version of B : phi^T*B

    CAROM::Vector reducedHn(numColumnRB, false);     // Reduced version of Hn
    reducedHn = 0;                                   // initialized to 0
    CAROM::Vector reducedHnm1(numColumnRB, false);   // Reduced version of Hn-1
    reducedHnm1 = 0;

    CAROM::Vector reducedBn(numColumnRB, false);     // Reduced version of Bn 
    reducedBn = 0;
    CAROM::Vector reducedBnm1(numColumnRB, false);   // Reduced version of Bn-1 
    reducedBnm1 = 0;

    CAROM::Vector reducedBdtn(numColumnRB, false);   // Reduced version of dBn/dt
    reducedBdtn = 0;

    CAROM::Vector reducedEx(numColumnRB, false);     // Reduced version of Exn
    reducedEx = 0;
    CAROM::Vector reducedEy(numColumnRB, false);     // Reduced version of Eyn
    reducedEy = 0;

    CAROM::Vector reducedExnm1(numColumnRB, false);  // Reduced version of Eyn-1
    reducedExnm1 = 0;
    CAROM::Vector reducedEynm1(numColumnRB, false);  // Reduced version of Eyn-1
    reducedEynm1 = 0;

    CAROM::Vector reducedJxnm1(numColumnRB, false);  // Reduced version of Jxn
    reducedJxnm1 = 0;
    CAROM::Vector reducedJynm1(numColumnRB, false);  // Reduced version of Jyn
    reducedJynm1 = 0;

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
        sout << "solution\n" << pmesh << Pn << "\nkeys j\n" << std::flush;

    }


    int vis_steps = 100;  // Number of iterations between two visualizations
    
    // The boundary condition is projected into the reduced space 
    CAROM::Vector B_carom_bdr(B.GetData(), B.Size(), true, true);  
    std::unique_ptr<CAROM::Vector> B_bdr_reduced = spatialbasis->transposeMult(B_carom_bdr); // Will contain the real dbc
    std::unique_ptr<CAROM::Vector> B_bdr_inv_r = spatialbasis->transposeMult(B_carom_bdr);   // Will contain only 1/(2*pi*r)

    // LinearForm to compute the flux by doing lf_flux^T . Bn = integral of Bn = flux_n
    ConstantCoefficient ones(1.);
    LinearForm lf_flux(fespace);
    lf_flux.AddDomainIntegrator(new DomainLFIntegrator(ones));
    lf_flux.Assemble();

    // LinearForm to compute the integral of the power density
    LinearForm lf_p(fespace);
    lf_p.AddDomainIntegrator(new DomainLFIntegrator(r_coeff));  // Adding r coefficient as the integral is on the volume in cylindrical coordinates
    lf_p.Assemble();

    // To be able to compute the power losses directly on the reduced space
    // It is possible to express the Power losses as J^T . Mp . E where Mp is
    // a diagonal matrix such as diag(Mp) = lf_p defined just above.
    // This method allows us to project the matrix in the reduced space and then
    // to directly use the reduced version of J and E.
    SparseMatrix Matrix_power(lf_p);
    CAROM::Matrix Reduced_Matrix_power_X(numColumnRB, numColumnRB, false);
    CAROM::Matrix Reduced_Matrix_power_Y(numColumnRB, numColumnRB, false);

    // Computing the projection of Mp defined just above in the reduced space for each component
    ComputeCtAB_Sparse(Matrix_power, *spatialbasisEx, *spatialbasisEx, Reduced_Matrix_power_X);
    ComputeCtAB_Sparse(Matrix_power, *spatialbasisEy, *spatialbasisEy, Reduced_Matrix_power_Y);


    // **** Computing the lf_flux in the reduced space ****
    // reduced_lf_flux = SpatialBasis^T . lf_flux 
    CAROM::Vector reduced_lf_flux(numColumnRB, false);
    CAROM::Vector lf_flux_carom(lf_flux.GetData(), lf_flux.Size(), true, false);
    spatialbasis->transposeMult(lf_flux_carom, reduced_lf_flux);

    auto t0 = std::chrono::high_resolution_clock::now(); // Time at the begining of the simulation
    for (int step = 0; step < num_steps; step++)
    {   
        iter++;

        t += Ts;  // Update the time

        // As the problem is linear, we can directly rescale the boundary conditions in the reduced basis
        *B_bdr_reduced = *B_bdr_inv_r;  // dbc : 1/(2*pi*r) in the reduced space
        *B_bdr_reduced *= NI_func(t);   // dbc : NI/(2*pi*r) in the reduced space
    
        // reducedRHS = reducedR2 * reduced_hnm1 + reducedM * reducedDBdT
        std::unique_ptr<CAROM::Vector> reduced_R2Hnm1 = ReducedR2.mult(reducedHnm1);
        std::unique_ptr<CAROM::Vector> reduced_MdBdt = ReducedM.mult(reducedBdtn);
        *reducedRHS = *reduced_R2Hnm1;
        *reducedRHS += *reduced_MdBdt;

        // Adding dbc effects in the RHS
        *reducedRHS += *B_bdr_reduced;

        // Getting the solution ReducedHn = reducedA^(-1) * reducedRHS
        invReducedA.mult(*reducedRHS, reducedHn);


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


        std::unique_ptr<CAROM::Vector> reduced_Jx = ReducedCx.mult(reducedHn);
        std::unique_ptr<CAROM::Vector> reduced_Jy = ReducedCy.mult(reducedHn);

        reducedExnm1 *= A3;
        reducedJxnm1 *= A2;
        *reduced_Jx *= A1;
        reducedEx = *reduced_Jx;
        reducedEx += reducedJxnm1;
        reducedEx += reducedExnm1;

        reducedEynm1 *= A3;
        reducedJynm1 *= A2;
        *reduced_Jy *= A1;
        reducedEy = *reduced_Jy;
        reducedEy += reducedJynm1;
        reducedEy += reducedEynm1;

        *reduced_Jx *= 1/A1;
        *reduced_Jy *= 1/A1;
        
        reducedJxnm1 = *reduced_Jx;
        reducedJynm1 = *reduced_Jy;
        reducedExnm1 = reducedEx;
        reducedEynm1 = reducedEy;


        // ******* Computing the power lossses *******
        // P = reduced_J^T . reduced_Mp . reduced_E
        CAROM::Vector MEx(numColumnRB, false);
        Reduced_Matrix_power_X.mult(reducedEx, MEx);
        real_t Px = MEx.inner_product(*reduced_Jx);

        CAROM::Vector MEy(numColumnRB, false);
        Reduced_Matrix_power_Y.mult(reducedEy, MEy);
        real_t Py = MEy.inner_product(*reduced_Jy);

        real_t vol = M_PI * height * (Rout*Rout - Ri*Ri);
        real_t P_eddy_tot = 2*M_PI*(Px + Py);
        real_t P_eddy = P_eddy_tot/vol;

        real_t flux= reduced_lf_flux.inner_product(reducedBn);
        data_file << t << ";" << P_eddy << ";" << flux << ";" << NI_func(t) << std::endl;


        if (step % vis_steps == 0)
        {       
            std::cout << "Time step: " << step << ", Time: " << t << std::endl;
            if (visualization)
            {
                spatialbasis->mult(reducedHn, X_carom);
                spatialbasis->mult(reducedBn, Bn_carom);
                spatialbasisEx->mult(*reduced_Jx, jx_carom);
                spatialbasisEy->mult(*reduced_Jy, jy_carom);
                spatialbasisEx->mult(reducedEx, Ex_carom);
                spatialbasisEy->mult(reducedEy, Ey_carom);

                for (int i = 0; i<Hn.Size(); i++){       
                    Hn(i) = X_carom(i);                  // Getting Hn in a mfem::Vector
                    Bn(i) = Bn_carom(i);
                    Pn(i) = jx_carom(i) * Ex_carom(i) + jy_carom(i) * Ey_carom(i);
                    Jn(i) = jx_carom(i);
                    Jn(i + Hn.Size()) = jy_carom(i);
                    // Jnx(i) = jx_carom(i); // No need as data are shared 
                    // Jny(i) = jy_carom(i);
                    // Enx(i) = Ex_carom(i); // No need as data are shared 
                    // Eny(i) = Ey_carom(i);
                    En(i) = Ex_carom(i);
                    En(i + Hn.Size()) = Ey_carom(i);
                }

                sout.precision(8);
                sout << "solution\n"
                     << pmesh << Hn
                     << "window_title 'Champ H'"
                     << std::flush;
                std::cout << "Press to continue" << std::endl;
                std::cin.get();
            } 
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();     // Capture time after the simulation
    std::chrono::duration<double> dt = t1 - t0;              // Duration of the simulation
    std::cout << "Elapsed time for solving ROM: " <<  dt.count() << " s\n";
    std::cout << "ROM Simulation has ended succesfully." << std::endl;
    delete fespace_E;
    delete fespace;
    delete fec;
    delete Cx;
    delete Cy;
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