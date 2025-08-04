#include "solver.hpp"
#include <iostream>
#include <cmath>


using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif



// Boundary condition parameters



PowerLossCoefficient::PowerLossCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> rho_eq_, CurlCustomCoefficient &J_r_, CurlCustomCoefficient &J_i_)
    : fespace(fespace_),
      J_r(J_r_), J_r_vect(fespace_->GetMesh()->SpaceDimension()),
      J_i(J_i_), J_i_vect(fespace_->GetMesh()->SpaceDimension()),
      rho_eq(rho_eq_) {}

real_t PowerLossCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    J_r.Eval(J_r_vect, T, ip);    // Get from J_r (Coefficient) the value at the point ip in J_r_vect
    J_i.Eval(J_i_vect, T, ip);    // same for imaginary part
    return  (J_r_vect * J_r_vect + J_i_vect * J_i_vect) * rho_eq.real();  // Re(rho) * J² * r (Cylindrical coordinates)
}


PowerLossMagCoefficient::PowerLossMagCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> mu_eq_, real_t omega_,  GridFunctionCoefficient &H_r_, GridFunctionCoefficient &H_i_)
    : fespace(fespace_),
      H_r(H_r_), H_i(H_i_),
      mu_eq(mu_eq_), omega(omega_) {}

real_t PowerLossMagCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{

    std::complex<real_t> h(H_r.Eval(T, ip),H_i.Eval(T, ip));
    return -(omega * std::imag(mu_eq) * pow(std::abs(h), 2));  
}



// Function for computing the power loss (by eddy currents)
void GetPowerLoss(Mesh *mesh, real_t f, real_t Bpeak, real_t &NI, Ferrite ferrite, real_t &P_loss_eddy, real_t &P_loss_mag, std::complex<real_t> &flux, const bool visualization) {

    // ---------------- Material Properties --------------
    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;
    real_t mu = ferrite.mu;

    real_t fc_mu = ferrite.fc;
    real_t tau = 1./(2.*M_PI*fc_mu);    // cutoff frequency of µ_eq = µ0µr/(1 + tau*j*w) model

    
    // --------------- Complex values coefficients --------
    real_t omega = 2*M_PI*f;
    const std::complex<real_t> j(0., 1.);  // imaginary number such that j² = -1 

    const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps); 
    const std::complex<real_t> mu_eq = mu / ((real_t)1. + tau * j * omega);
    const std::complex<real_t> j_omega_mu = j*omega*mu_eq;

    std::cout << "|rho| = " << std::abs(rho_eq) << std::endl;
    // The value of the source NI is initialized to 1. As the problem is linear, we can later rescale the solution.
    // We are doing such a thing because the actual source is the magnetic flux, but we don't know the corresponding current
    NI = 1;

    // Boundary conditions function
    auto bdr_func = [NI](const Vector &x)
    {
        real_t r = x[0];
        return NI /(2 * M_PI * r);
    };

    // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];  
    };

    // Function for the "correcting" factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };


    int order = 1;                      // Order of shape functions
    int dim = mesh->Dimension();        // Mesh dimension (2D)

    // Working with scalar values as it is an axisymetric problem
    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    // Since we are in the harmonic steady-state regime, we work with complex values
    ComplexGridFunction h(fespace);

    // Definition of the Boundary conditions  :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // All borders have Dirichlet boundary conditions (dbc), thus all border attributes are set to 1
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // Getting the list of degrees of freedom with dbc

    // FunctionCoefficient for the boundary conditions
    FunctionCoefficient bdr_coeff(bdr_func);  // Boundary condition of NI/(2*pi*r) on the border
    real_t zero = 0.;
    ConstantCoefficient zero_coeff(zero);     // Used for the complexe value

    // Projecting the boundary conditions on the border elements of the GridFunction h
    h.ProjectBdrCoefficient(bdr_coeff, zero_coeff, dir_bdr); // NI/(2*pi*r) for real part and 0 for imaginary part
    
    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Hermitian convention for the hermitian product
    b.Assemble();  // LinearForm equals to 0 as there is no source term

    // Coefficients for the definition of the bilinear form

    FunctionCoefficient r_coeff(r_coeff_func);   // Coefficient r for the integral in cylindrical coordinates

    // rho_eq Coefficient separated in real and imaginary part 
    ConstantCoefficient rho_real(rho_eq.real()); // for the -div(rho grad(h)) terme
    ConstantCoefficient rho_imag(rho_eq.imag());
    
    // **** Coefficient for the integration in cylindrical coordinates 
    ProductCoefficient r_rho_real(rho_real, r_coeff);   // r * rho_r
    ProductCoefficient r_rho_imag(rho_imag, r_coeff);   // r * rho_i
    // *****

    // j*omega*mu coefficient separated in real and imaginary part
    ConstantCoefficient j_omega_mu_real(j_omega_mu.real());  // jw * mu_r
    ConstantCoefficient j_omega_mu_imag(j_omega_mu.imag());  // jw * mu_i
    
    // **** Coefficient j*omega*mu_eq*r for the integration in cylindrical coordinates 
    ProductCoefficient r_j_omega_mu_real(j_omega_mu_real, r_coeff);  // r * jw * mu_r
    ProductCoefficient r_j_omega_mu_imag(j_omega_mu_imag, r_coeff);  // r * jw * mu_i
    // ****

    // rho_eq/r² term for correcting the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);                  // 1/r
    ProductCoefficient rho_inv_square_coeff_r(inv_r_square_coeff, r_rho_real);  // r * rho_r/r²
    ProductCoefficient rho_inv_square_coeff_i(inv_r_square_coeff, r_rho_imag);  // r * rho_i/r²
    

    // As we are using complex values, it is not a bilinear form but a sesquilinear form
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN); // Hermitian convention
    a->AddDomainIntegrator(new DiffusionIntegrator(r_rho_real), new DiffusionIntegrator(r_rho_imag)); //  -div(rho_eq grad(H_theta))
    a->AddDomainIntegrator(new MassIntegrator(r_j_omega_mu_real), new MassIntegrator(r_j_omega_mu_imag));             //  jw mu_eq H

    // rho_eq/r² term for "correcting" the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    a->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r), new MassIntegrator(rho_inv_square_coeff_i));  
    a->Assemble();
    
    // Defining a preconditioner using the real part
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(r_rho_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(r_j_omega_mu_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r));
    pcOp->Assemble();

    OperatorHandle A;
    Vector B, H; 
    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B); 

    // Setting a preconditioner 
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

    pc_r = new GSSmoother(*PCOp.As<SparseMatrix>());  // Using a Gauss-Seidel smoother

    real_t s = -1.0;
    pc_i = new ScaledOperator(pc_r, s);  // Same preconditioner for omaginary part, only multiplied by -1 

    BDP.SetDiagonalBlock(0, pc_r);
    BDP.SetDiagonalBlock(1, pc_i);
    BDP.owns_blocks = 1;

    // --------- Solver ----------
    GMRESSolver gmres;
    gmres.SetPreconditioner(BDP);  // Applying the preconditioner
    gmres.SetOperator(*A.Ptr()); 
    gmres.SetRelTol(1e-12);
    gmres.SetMaxIter(2000);
    gmres.SetPrintLevel(0);
    gmres.Mult(B, H);   // Solve the system

    // 
    a->RecoverFEMSolution(H,b,h);  // Get system solution in h ComplexGridFunction
   
    GridFunction h_r = h.real();  // Getting real and imaginary part
    GridFunction h_i = h.imag();

    // GridFunction used for computing integral on the domain (for the magnetic flux)
    GridFunction ones(fespace);
    ones = 1; // All values set to 1

    ConstantCoefficient one_coeff(1.);

    // Coefficients for computing the integral of B (magnetic flux)
    GridFunctionCoefficient H_r_coeff(&h_r);
    GridFunctionCoefficient H_i_coeff(&h_i);


    // LinearForm representing the integral of H
    LinearForm lf_flux(fespace);
    lf_flux.AddDomainIntegrator(new DomainLFIntegrator(one_coeff));
    lf_flux.Assemble();
    real_t flux_r = lf_flux(h_r);
    real_t flux_i = lf_flux(h_i);

    flux = flux_r + j * flux_i;  // Flux of H
    flux *= mu_eq;               // Flux of B

    real_t flux_abs = std::abs(flux);

    GridFunction B_average(fespace);    // Constant GridFunction containing the average targeted value of B
    B_average = Bpeak/sqrt(2);          // Working with RMS values

    real_t targeted_flux = lf_flux(B_average);

    real_t coeff_correction;
    coeff_correction = targeted_flux / flux_abs;

    // Adjusting h values such that the magnetic flux is equal to the desired value
    h_r *= coeff_correction;
    h_i *= coeff_correction;

    // Computing the coresponding NI 
    NI *= coeff_correction; // RMS value
    flux *= coeff_correction;


    // std::cout << "Corrective factor : " << coeff_correction << std::endl;
    // std::cout << "NI (peak value) = " << NI*sqrt(2) << std::endl;



    // ************* Eddy losses *****************
    // Using L2 finite elements as the curl of a H1-space is discountinuous
    FiniteElementCollection *fec_l2 = new L2_FECollection(order, dim);

    FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec_l2, dim);    // Space for Current Density / Electrical field
    FiniteElementSpace *fespace_P = new FiniteElementSpace(mesh, fec_l2);         // Field for the power losses 

    // Computing curl(h) as a coefficient (curl(H) = J)
    CurlCustomCoefficient curl_r_coeff(&h_r);
    CurlCustomCoefficient curl_i_coeff(&h_i);

    GridFunction ones_L2(fespace_P);  // GridFunction For the volume
    ones_L2 = 1;

    // using the Coefficient curl(h), it is possible to compute the power loss density (Re(rho_eq) * (J_r² + J_i²))
    PowerLossCoefficient Power_loss_coef(fespace_P, rho_eq, curl_r_coeff, curl_i_coeff);
    GridFunction P_eddy(fespace_P);
    P_eddy.ProjectCoefficient(Power_loss_coef);

    // Integrating the power loss density
    LinearForm lf_vol(fespace_P);
    lf_vol.AddDomainIntegrator(new DomainLFIntegrator(r_coeff));
    lf_vol.Assemble();

    real_t vol = lf_vol(ones_L2);

    P_loss_eddy = lf_vol(P_eddy)/vol; // Power in Watts / m^3
    

    // ************* Mag losses *****************
    PowerLossMagCoefficient Power_loss_mag_coeff(fespace_P, mu_eq, omega, H_r_coeff, H_i_coeff);
    GridFunction P_mag(fespace_P);
    P_mag.ProjectCoefficient(Power_loss_mag_coeff);

    P_loss_mag = lf_vol(P_mag)/vol; // Power in Watts / m^3
          
    // std::cout << "f = " << f << " Hz, " << "Ploss : " << P_loss_eddy + P_loss_mag << " W/m^3"<< std::endl;
    // std::cout << "P_eddy = " << P_loss_eddy << " W/m^3"<< std::endl;
    // std::cout << "P_mag = " << P_loss_mag << " W/m^3"<< std::endl;
    if (visualization){

        h_r *= sqrt(2); // Plotting peak values
        h_i *= sqrt(2);

        //GridFunctions for J (only used for plot)
        GridFunction J_r(fespace_E);
        GridFunction J_i(fespace_E);

        //Projecting Coefficient values on the GridFunctions
        J_r.ProjectCoefficient(curl_r_coeff);
        J_i.ProjectCoefficient(curl_i_coeff);

        // Visualisation GLVis
        char vishost[] = "localhost";
        int  visport   = 19916;
        socketstream sol_sock_r(vishost, visport);
        socketstream sol_sock_i(vishost, visport);
        socketstream sol_sock_j_r(vishost, visport);
        socketstream sol_sock_j_i(vishost, visport);

        sol_sock_r.precision(8);
        sol_sock_i.precision(8);
        sol_sock_j_r.precision(8);
        sol_sock_j_i.precision(8);

        sol_sock_r << "solution\n" << *mesh << h_r
                    << "window_title 'Solution: Real Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        sol_sock_i << "solution\n" << *mesh << h_i
                    << "window_title 'Solution: Imaginary Part'" 
                    << "pause\n" << "keys c\n" << std::flush;

        sol_sock_j_r << "solution\n" << *mesh << J_r
                    << "window_title 'Solution: J Real Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        sol_sock_j_i << "solution\n" << *mesh << J_i
                    << "window_title 'Solution: J Imaginary Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        
    }
    // Free memory
    delete fec;
    delete fespace;
    delete a;
    delete pcOp;
    delete fespace_E;
}