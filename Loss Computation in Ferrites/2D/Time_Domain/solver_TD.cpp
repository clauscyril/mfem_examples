#include "solver_TD.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

void TD_sim(Mesh *mesh, real_t I_rms, real_t f, real_t Ts, int num_steps, Ferrite ferrite, bool visualization){

    // geometric parameters
    real_t Ri = 9.6e-3/2.0;
    real_t height = 7.59e-3;
    real_t w = 5.3e-3;
    real_t Rout = Ri + w;
    real_t Rm = (Rout - Ri) / 2;


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

    // // Boundary conditions function
    auto bdr_func = [&](const Vector &x) -> real_t
    {
        real_t r = x(0); // Accès au premier élément du vecteur x
        real_t NI = I_rms * std::sqrt(2) * std::sin(2 * M_PI * f * t);
        return NI / (2 * M_PI * r);
    };

    // // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];
    };


    // Function for the correcting factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };
    
    int order = 2;                      // elements order : 2
    int dim = mesh->Dimension();        // Mesh dimension : 2

    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec, dim);

    GridFunction Hn(fespace);     // H^n
    GridFunction Hnm1(fespace);   // H^{n-1}
    GridFunction db_dt(fespace);
    GridFunction db_dtm1(fespace);

    GridFunction Bn(fespace);
    GridFunction Bnm1(fespace);

    GridFunction En(fespace_E);
    GridFunction Enm1(fespace_E);

    GridFunction Jn(fespace_E);
    GridFunction Jnm1(fespace_E);
    // CurlCustomCoefficient *curl_H_coeff = new CurlCustomCoefficient(&Hn);

    GridFunction Pn(fespace);

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
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // All the borders have boundary conditions
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);

    // FunctionCoefficient bdr_coeff([&](const Vector &x) { return bdr_func(x, t); });
    FunctionCoefficient bdr_coeff(bdr_func);
    Hnm1.ProjectBdrCoefficient(bdr_coeff, dir_bdr);

    // **** Coefficients for bilinear forms
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


    BilinearForm r1(fespace); 
    r1.AddDomainIntegrator(new DiffusionIntegrator(A1_r));
    r1.AddDomainIntegrator(new MassIntegrator(B1_r));
    r1.AddDomainIntegrator(new MassIntegrator(A1_r_inv));
    r1.Assemble();

    BilinearForm r2(fespace); 
    r2.AddDomainIntegrator(new DiffusionIntegrator(A2_r));
    r2.AddDomainIntegrator(new MassIntegrator(B2_r));
    r2.AddDomainIntegrator(new MassIntegrator(A2_r_inv));
    r2.Assemble();


    // Initialisation à 0
    // GridFunctionCoefficient db_dt_coeff(&db_dt); 
    // LinearForm *m = new LinearForm(fespace);
    // m->AddDomainIntegrator(new DomainLFIntegrator(k3_r));
    // m->Assemble();
    BilinearForm m(fespace);
    m.AddDomainIntegrator(new MassIntegrator(k_r));
    m.Assemble();


    r1.Finalize();
    r2.Finalize();
    m.Finalize();


    SparseMatrix &A = r1.SpMat();
    Vector rhs(fespace->GetTrueVSize());
    Vector R2Hnm1(fespace->GetTrueVSize());
    Vector MdBdt(fespace->GetTrueVSize());

    Vector x(fespace->GetTrueVSize()); // Pour stocker solution complète

    CGSolver solver;
    solver.SetOperator(A);
    solver.SetRelTol(1e-12);
    solver.SetMaxIter(5000);
    solver.SetPrintLevel(0);

    socketstream sout;
    socketstream sout_e;
    socketstream sout_j;
    socketstream sout_b;

    char vishost[] = "localhost";
    int visport = 19916;

    if (visualization)
    {
        sout.open(vishost, visport);
        sout_e.open(vishost, visport);
        sout_j.open(vishost, visport);
        sout_b.open(vishost, visport);

        sout.precision(8);
        sout << "solution\n" << *mesh << Hn << "\nkeys j\n" << std::flush;
        sout_e.precision(8);
        sout_e << "solution\n" << *mesh << En << "\nkeys j\n" << std::flush;
        sout_j.precision(8);
        sout_j << "solution\n" << *mesh << Jn << "\nkeys j\n" << std::flush;
        sout_b.precision(8);
        sout_b << "solution\n" << *mesh << Bn << "\nkeys j\n" << std::flush;
    }

    std::string name = "./data/TD_" + std::to_string((int)f/1000) + ".csv";   // Path to csv file for python plot
    std::ofstream data_file(name);                          // ofstream for writing in the file
    data_file << "t;p_eddy;flux\n0;0;0\n";                    // Intialising the file with coluns names



    int vis_steps = 1000;
    // ***************  Time iterations *****************
    for (int step = 0; step < num_steps; step++)
    {
        t += Ts;
        bdr_coeff.SetTime(t);       // Update boundary conditions fonction's time parameters
        Hn.ProjectBdrCoefficient(bdr_coeff, dir_bdr); // Set new boundary conditions on H

        // Mult M * dBdt → MdBdt
        m.Mult(db_dt, MdBdt);

        // Mult R2 * Hnm1 → R2Hnm1
        r2.Mult(Hnm1, R2Hnm1);

        // RHS = R2Hnm1 + MdBdt
        rhs = R2Hnm1;
        rhs += MdBdt;
        
        SparseMatrix A_sys;
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

        GridFunction H_temp(fespace); // temporary GridFunction to avoid modifying Hn
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
        GridFunction one(fespace);
        one = 1;
        LinearForm lf_flux(fespace);
        lf_flux.AddDomainIntegrator(new DomainLFIntegrator(B_coeff));
        lf_flux.Assemble();

        real_t flux = lf_flux(one);

        // Computation of J and E
        CurlCustomCoefficient J_coeff(&Hn);
        Jn.ProjectCoefficient(J_coeff);  // Jn = curl(Hn)

        GridFunction Jn_temp(fespace_E);
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
        LinearForm lf(fespace);
        lf.AddDomainIntegrator(new DomainLFIntegrator(P_eddy_coeff));
        lf.Assemble();
        real_t P_eddy_tot = 2*M_PI*lf(one);
        real_t vol = M_PI * height * (Rout*Rout - Ri*Ri);
        real_t P_eddy = P_eddy_tot/vol;

        Pn.ProjectCoefficient(P_eddy_coeff);


        data_file << t << ";" << P_eddy << ";" << flux << std::endl;

        // Visualisation avec GLVis
        if (step % vis_steps == 0 )
        {
            std::cout << "Time step: " << step << ", Time: " << t << std::endl;

            if (visualization)
            {
                sout.precision(8);
                sout << "solution\n"
                     << *mesh << Hn
                     << "window_title 'Champ H'"
                     << std::flush;
                sout_e << "solution\n" << *mesh << En <<  "window_title 'Champ E'" << std::flush;
                sout_j << "solution\n" << *mesh << Jn <<  "window_title 'Courant J'" << std::flush;
                sout_b << "solution\n" << *mesh << Bn <<  "window_title 'Field B'" << std::flush;

                // std::cin.get(); // Décommenter pour pause manuelle
            }
        }
        
    }
    data_file.close();

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
    E.Eval(E_vect, T, ip);    // same
    return  (J_vect*E_vect) * x[0];  // Re(rho) * J² * r (Cylindrical coordinates)
}

Ferrite::Ferrite(const char* name_, real_t rho_, real_t sigma_, real_t eps_, real_t mu_r_)
    : name(name_),  rho(rho_), sigma(sigma_), eps(eps_), mu(4e-7*M_PI*mu_r_){}