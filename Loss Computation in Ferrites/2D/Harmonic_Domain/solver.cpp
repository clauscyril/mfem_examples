#include "solver.hpp"
#include <iostream>
#include <cmath>


using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

// geometric parameters
real_t Ri = 9.6e-3/2.0;
real_t height = 7.59e-3;
real_t w = 5.3e-3;
real_t Rout = Ri + w;
real_t Rm = (Rout - Ri) / 2;

// ****** Ferrite parameters ******

// Ferrite N30
real_t rho = 5.98e-2;
real_t sigma = 4.44e-1;
real_t eps = 2.48e-6;
real_t mu = 4300.0 * 4e-7 * M_PI;

// // Ferrite N87
// real_t rho = 4.24e-2;
// real_t sigma = 1.48e-1;
// real_t eps = 2.68e-6;
// real_t mu = 2200.0 *4e-7 * M_PI;

// *********************************

// Boundary condition parameters
real_t Bpeak = 10e-3;   // average value of the amplitude of the magnetic field through the surface (flux imposed)
real_t NI;              // Equivalent value of NI with Bpeak input (Needs to be rescaled in post processing)


PowerLossCoefficient::PowerLossCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> rho_eq_, CurlCustomCoefficient &J_r_, CurlCustomCoefficient &J_i_)
    : fespace(fespace_),
      J_r(J_r_), J_r_vect(fespace_->GetMesh()->SpaceDimension()),
      J_i(J_i_), J_i_vect(fespace_->GetMesh()->SpaceDimension()),
      rho_eq(rho_eq_) {}

real_t PowerLossCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    Vector x;               // Vector coordinates
    T.Transform(ip, x);     // Get the global coordinates in vector x from integration point's coordinates in the element referential
    J_r.Eval(J_r_vect, T, ip);    // Get from J_r (Coefficient) the value at the point ip in J_r_vect
    J_i.Eval(J_i_vect, T, ip);    // same
    return  (J_r_vect * J_r_vect + J_i_vect * J_i_vect) * rho_eq.real() * x[0];  // Re(rho) * J² * r (Cylindrical coordinates)
}


PowerLossMagCoefficient::PowerLossMagCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> mu_eq_, real_t omega_,  GridFunctionCoefficient &H_r_, GridFunctionCoefficient &H_i_)
    : fespace(fespace_),
      H_r(H_r_), H_i(H_i_),
      mu_eq(mu_eq_), omega(omega_) {}

real_t PowerLossMagCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    Vector x;               // Vector coordinates
    T.Transform(ip, x);     // Get the global coordinates in vector x from integration point's coordinates in the element referential
    std::complex<real_t> h(H_r.Eval(T, ip),H_i.Eval(T, ip));
    return -(omega * std::imag(mu_eq) * pow(std::abs(h), 2) * x[0]);  
}



// Function for computing the power loss (by eddy currents)
void GetPowerLoss(Mesh *mesh, real_t fc, real_t fc_mu, real_t &P_loss_eddy, real_t &P_loss_mag, std::complex<real_t> &flux, real_t &Imax, const bool visualization) {
    
    real_t omega = 2*M_PI*fc;           // Working frequency
    real_t tau = 1./(2.*M_PI*fc_mu);    // cutoff frequency of µ_eq = µ0µr/(1 + tau*j*w) model

    // Complexe parameters 
    const std::complex<real_t> j(0., 1.);  // j² = -1

    // Physical parameters model
    const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps); 
    const std::complex<real_t> mu_eq = mu / ((real_t)1. + tau * j * omega);

    const std::complex<real_t> j_omega_mu = j*omega*mu_eq;

    NI = 1.;

    // Boundary conditions function
    auto bdr_func = [](const Vector &x)
    {
        real_t r = x[0];
        return NI /(2 * M_PI * r); // IMPORTANT : NI must be defined in the code before calling this function
    };

    // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];  
    };

    // Function for the correcting factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };


    std::cout <<"NI : " <<  NI << ", |mu_eq| = " << std::abs(mu_eq) << std::endl;


    // working with 1D elements in a 2D space
    int order = 1;                      // elements order : 
    int dim = mesh->Dimension();        // Mesh dimension : 2
    std::cout << "test : " << dim << std::endl;

    // Working with scalar values as it is an axisymetric problem
    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    // Since we are in the harmonic steady-state regime, we work with complex values
    ComplexGridFunction h(fespace);

    // Definition of the Boundary conditions  :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // All the borders have boundary conditions
    // dir_bdr[1] = 1; 
    // dir_bdr[2] = 1; 
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);   // ess_tdof_list is a list of the physical regions containing Dirichlet Boundary condition 
                                                           
    // FunctionCoefficient for the boundary conditions
    FunctionCoefficient bdr_coeff(bdr_func);
    real_t zero = 0.;
    ConstantCoefficient zero_coeff(zero);

    // // Prjecting the boundary conditions on the GridFunction h
    h.ProjectBdrCoefficient(bdr_coeff, zero_coeff, dir_bdr); // NI/(2*pi*r) for real part and 0 for imaginary part
    
    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Hermitian convention for the hermitian product
    b.Assemble();  // LinearForm equals to 0 as there is no source term

    // Coefficients for the definition of the bilinear form

    FunctionCoefficient r_coeff(r_coeff_func);   // Coefficient r for the integral in cylindrical coordinates

    // rho_eq Coefficient separated in real and imaginary part 
    ConstantCoefficient rho_real(rho_eq.real()); // for the -div(rho grad(h)) terme
    ConstantCoefficient rho_imag(rho_eq.imag());
    
    // **** Coefficient for the integration in cylindrical coordinates 
    ProductCoefficient r_rho_real(rho_real, r_coeff);
    ProductCoefficient r_rho_imag(rho_imag, r_coeff);
    // *****

    // j*omega*mu coefficient separated in real and imaginary part
    ConstantCoefficient j_omega_mu_real(j_omega_mu.real()); 
    ConstantCoefficient j_omega_mu_imag(j_omega_mu.imag());
    
    // **** Coefficient j*omega*mu_eq*r for the integration in cylindrical coordinates 
    ProductCoefficient r_j_omega_mu_real(j_omega_mu_real, r_coeff);
    ProductCoefficient r_j_omega_mu_imag(j_omega_mu_imag, r_coeff);
    // ****

    // rho_eq/r² term for correcting the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);
    ProductCoefficient rho_inv_square_coeff_r(inv_r_square_coeff, r_rho_real);
    ProductCoefficient rho_inv_square_coeff_i(inv_r_square_coeff, r_rho_imag);
    

    // ComplexBilinearForm (SesquilinearForm)
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    a->AddDomainIntegrator(new DiffusionIntegrator(r_rho_real), new DiffusionIntegrator(r_rho_imag)); //  -div(rho_eq grad(H_theta))
    a->AddDomainIntegrator(new MassIntegrator(r_j_omega_mu_real), new MassIntegrator(r_j_omega_mu_imag));             //  jw mu_eq H

    // rho_eq/r² term for correcting the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    a->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r), new MassIntegrator(rho_inv_square_coeff_i));  

    a->Assemble();
    
    // Defining a preconditioner using the real part
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(r_rho_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(r_j_omega_mu_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r));

    pcOp->Assemble();

    OperatorHandle A, Pc;
    Vector B, H; 
    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B); 

    // Using the preconditioner `pc0p` in the solver
    pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    GSSmoother M((SparseMatrix&)(*Pc));    // Using a Gauss-Seidel smoother

    Array<int> blockOffsets;
    blockOffsets.SetSize(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = A->Height() / 2;
    blockOffsets[2] = A->Height() / 2;
    blockOffsets.PartialSum();

    BlockDiagonalPreconditioner BDP(blockOffsets);

    Operator * pc_r = NULL;
    Operator * pc_i = NULL;

    OperatorHandle PCOp;
    
    pcOp->SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
    pcOp->FormSystemMatrix(ess_tdof_list, PCOp);


    // // Il y a plusieurs types de préconditionneurs disponibles :
    // pc_r = new DSmoother(*PCOp.As<SparseMatrix>());

    pc_r = new GSSmoother(*PCOp.As<SparseMatrix>());

    // // Si on utilise les calculs parallèles (ParMesh, etc...), on peut utiliser celui-ci :
    // HypreBoomer   AMG *amg = new HypreBoomerAMG();
    // amg->SetOperator(*PCOp.As<HypreParMatrix>());
    // pc_r = amg;

    real_t s = -1.0;
    pc_i = new ScaledOperator(pc_r, s);  // On utilise le même preconditionneur pour la partie imaginaire, mais multiplié par -1

    BDP.SetDiagonalBlock(0, pc_r);
    BDP.SetDiagonalBlock(1, pc_i);

    BDP.owns_blocks = 1;

    GMRESSolver gmres;
    gmres.SetPreconditioner(BDP);
    gmres.SetOperator(*A.Ptr());
    gmres.SetRelTol(1e-12);
    gmres.SetMaxIter(100000);
    gmres.SetPrintLevel(0);
    gmres.Mult(B, H);

    
    // UMFPackSolver solver;
    // solver.SetOperator(*A.Ptr());
    // solver.Mult(B, H);

    a->RecoverFEMSolution(H,b,h);  // Get system solution in h GridFunction
   
    GridFunction h_r = h.real();  // Getting real and imaginary part
    GridFunction h_i = h.imag();

    // GridFunction used for computing integral on the domain (for the magnetic flux)
    GridFunction ones(fespace);
    ones = 1; // All values set to 1

    // GridFunctions for the magnetic field B 
    GridFunction B_r(fespace);
    GridFunction B_i(fespace);

    // Computing B with B = µH
    // B_r = h_r;
    // B_r *= mu_eq.real();

    // B_i = h_i;
    // B_i *= mu_eq.imag();

    // // Coefficients for computing the integral of B (magnetic flux)
    // GridFunctionCoefficient B_r_coeff(&B_r);
    // GridFunctionCoefficient B_i_coeff(&B_i);

        // Coefficients for computing the integral of B (magnetic flux)
    GridFunctionCoefficient H_r_coeff(&h_r);
    GridFunctionCoefficient H_i_coeff(&h_i);

    // LinearForm representing the integral of B
    LinearForm lf_r_flux(fespace);
    lf_r_flux.AddDomainIntegrator(new DomainLFIntegrator(H_r_coeff));
    lf_r_flux.Assemble();
    real_t flux_r = lf_r_flux(ones);

    // The same for the imaginary part
    LinearForm lf_i_flux(fespace);
    lf_i_flux.AddDomainIntegrator(new DomainLFIntegrator(H_i_coeff));
    lf_i_flux.Assemble();
    real_t flux_i = lf_i_flux(ones);

    flux = flux_r + j * flux_i;
    flux *= mu_eq;
    real_t flux_abs = std::abs(flux);

    // Computing the corrective term to have the desired magnetic flux
    real_t B_rms;
    B_rms = flux_abs / w / height;
    std::cout << "Brms pre-processing = " << B_rms << ", B_peak pre-processing = " << B_rms*sqrt(2) << std::endl;
    // Printing the magnetic flux
    std::cout << "Flux magnétique = " << flux_r << " + " << flux_i << "j" << std::endl;
    std::cout << "Flux magnétique (module) = " << flux_abs << std::endl;

    real_t coeff_correction;
    if (Imax == 0.)  // Imax = 0 si on impose avec Bpeak
    {
        coeff_correction = Bpeak/sqrt(2) / B_rms;
        Imax = sqrt(2) * NI * coeff_correction;
    }
    else
    {
        coeff_correction = Imax/(Bpeak/mu *M_PI * Rm);
    }

    // Adjusting h values such that the magnetic flux is equal to the desired value
    h_r *= coeff_correction;
    h_i *= coeff_correction;

    // Computing the coresponding NI 
    NI *= coeff_correction; // RMS value
    Imax = NI * sqrt(2);
    flux *= coeff_correction;
    B_rms *= coeff_correction;

    std::cout << "Brms post-processing = " << B_rms << ", B_peak post-processing = " << B_rms*sqrt(2) << std::endl;

    // Printing the corrective term
    std::cout << "Terme correctif : " << coeff_correction << std::endl;
    std::cout << "NI (amplitude) = " << NI*sqrt(2) << std::endl;
    std::cout << "rho_eq = " << rho_eq.real() << " + j" << rho_eq.imag() << std::endl;
    std::cout << "mu_eq  = " << mu_eq.real() << " + j" << mu_eq.imag() << std::endl;


    // **********************************Power Losses ************************
    real_t vol = M_PI * height * (pow(Rout,2) - pow(Ri,2));

    // ************* Eddy losses *****************
    // Defining a new space for the current density J
    FiniteElementCollection *fec_E = new ND_FECollection(order, dim);
    FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec_E);
    // FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec, dim);

    // Computing curl(h) as a coefficient (curl(H) = J)
    CurlCustomCoefficient curl_r_coeff(&h_r);
    CurlCustomCoefficient curl_i_coeff(&h_i);

    //GridFunctions for J
    GridFunction J_r(fespace_E);
    GridFunction J_i(fespace_E);

    //Projecting Coefficient values on the GridFunctions
    J_r.ProjectCoefficient(curl_r_coeff);
    J_i.ProjectCoefficient(curl_i_coeff);

    // using the Coefficient curl(h), it is possible to compute the power loss density (Re(rho_eq) * (J_r² + J_i²) * r)
    PowerLossCoefficient Power_loss_coef(fespace_E, rho_eq, curl_r_coeff, curl_i_coeff);

    // Integrating the power loss density (remembering that we are in cylindrical coordinates)
    LinearForm lf_eddy(fespace);
    lf_eddy.AddDomainIntegrator(new DomainLFIntegrator(Power_loss_coef));
    lf_eddy.Assemble();
    real_t P_loss_eddy_tot = 2*M_PI*lf_eddy(ones); // Power in Watts

    P_loss_eddy = P_loss_eddy_tot/vol;  // Average Power in W/m^3

    std::cout << "Ploss : " << P_loss_eddy_tot << std::endl;
    std::cout << "Ploss(W/m^3) : " << P_loss_eddy << std::endl;

    // ************* Mag losses *****************
    // GridFunctionCoefficient H_r_coeff(&h_r);
    // GridFunctionCoefficient H_i_coeff(&h_i);

    PowerLossMagCoefficient Power_loss_mag_coeff(fespace_E, mu_eq, omega, H_r_coeff, H_i_coeff);

    LinearForm lf_mag(fespace);
    lf_mag.AddDomainIntegrator(new DomainLFIntegrator(Power_loss_mag_coeff));
    lf_mag.Assemble();
    real_t P_loss_mag_tot = 2*M_PI*lf_mag(ones); // Power in Watts
    P_loss_mag = P_loss_mag_tot/vol;




    
    h_r *= sqrt(2);
    h_i *= sqrt(2);
          
    if (visualization){
        // Visualisation GLVis
        char vishost[] = "localhost";
        int  visport   = 19916;
        socketstream sol_sock_r(vishost, visport);
        socketstream sol_sock_i(vishost, visport);
        // socketstream sol_sock_mod(vishost, visport);
        // socketstream sol_sock_phase(vishost, visport);
        sol_sock_r.precision(8);
        sol_sock_i.precision(8);
        // sol_sock_mod.precision(8);
        // sol_sock_phase.precision(8);
        sol_sock_r << "solution\n" << *mesh << h_r
                    << "window_title 'Solution: Real Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        sol_sock_i << "solution\n" << *mesh << h_i
                    << "window_title 'Solution: Imaginary Part'" 
                    << "pause\n" << "keys c\n" << std::flush;


        socketstream sol_sock_j_r(vishost, visport);
        socketstream sol_sock_j_i(vishost, visport);
        sol_sock_j_r.precision(8);
        sol_sock_j_i.precision(8);
        sol_sock_j_r << "solution\n" << *mesh << J_r
                    << "window_title 'Solution: J Real Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        sol_sock_j_i << "solution\n" << *mesh << J_i
                    << "window_title 'Solution: J Imaginary Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        
    }
    // ParaViewDataCollection paraview_dc("torre_solution", mesh);
    // paraview_dc.SetPrefixPath("paraview_output"); // Optionnel : dossier de sortie
    // paraview_dc.RegisterField("H_real", &h_r);
    // paraview_dc.RegisterField("H_imag", &h_i);
    // paraview_dc.SetLevelsOfDetail(order);
    // paraview_dc.SetHighOrderOutput(true);
    // paraview_dc.SetCycle(0);
    // paraview_dc.SetTime(0.0);
    // paraview_dc.Save();

    delete fec;
    delete fespace;
    delete a;
    delete pcOp;

    delete fespace_E;
}