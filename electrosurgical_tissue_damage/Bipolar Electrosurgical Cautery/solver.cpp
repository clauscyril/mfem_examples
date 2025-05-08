#include "solver.hpp"
#include <iostream>
#include <cmath>


using namespace mfem;


void compute_V(Mesh &mesh, FiniteElementSpace *fespace, Array<int> &ess_tdof_list, GridFunction &v){
    int affichage = 1;
    std::cout << "Computing q using fem" << std::endl;


    LinearForm b(fespace);
    // b = 0;
    b.Assemble();

    BilinearForm a(fespace);
    ConstantCoefficient sigma(0.2);
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    
    OperatorPtr A;
    Vector B, V;

    a.FormLinearSystem(ess_tdof_list, v, b, A, V, B);
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, V, 1, 1000, 1e-12, 0.f);

    // On récupère les solutions 
    a.RecoverFEMSolution(B, b, v);




    GradientGridFunctionCoefficient grad_v_coeff(&v);

    FiniteElementCollection *fec_grad = new RT_FECollection(1, 2);  // Raviart-Thomas (RT) pour le gradient
    FiniteElementSpace *fespace_grad = new FiniteElementSpace(&mesh, fec_grad);
    
    GridFunction grad_v(fespace_grad);

    grad_v.ProjectCoefficient(grad_v_coeff);
    GridFunction E(fespace_grad);
    E -= grad_v;
    E *= 0.258;









    // QCoefficient q_coeff(v, 0.22+0.56);
    // // ProductCoefficient product(&grad_v, &grad_v);

    // q.ProjectCoefficient(q_coeff);

    // Visualisation GLVis
    if (affichage){
        char vishost[] = "localhost";
        int  visport   = 19916;
        socketstream sol_sock_v(vishost, visport);
        sol_sock_v.precision(8);
        sol_sock_v << "solution\n" << mesh << v
                    << "window_title 'Solution: Potentiel'" 
                    << "pause\n" << "keys c\n" << std::flush;

        socketstream sol_sock_q(vishost, visport);
        sol_sock_q.precision(8);
        sol_sock_q << "solution\n" << mesh << E
                    << "window_title 'Solution: E'" 
                    << "pause\n" << "keys c\n" << std::flush;
    }
        
}


QCoefficient::QCoefficient(const GridFunction &gf_, real_t period_, real_t t_on_)
    : gf(gf_), grad_coeff(&gf_), 
      grad(gf_.FESpace()->GetMesh()->SpaceDimension()),
      time(0.0), period(period_), t_on(t_on_) {}

void QCoefficient::SetTime(real_t t) {
    time = t;
}

real_t QCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    grad_coeff.Eval(grad, T, ip);
    real_t norml2 = grad.Norml2();
    norml2 *= norml2;

    float phase = fmod(time, period);
    if (phase < t_on) 
        return 0.258 * norml2;
    return 0;
}

DG_Solver::DG_Solver(SparseMatrix &M_, SparseMatrix &K_, const FiniteElementSpace &fes)
    : M(M_), K(K_),
      prec(fes.GetTypicalFE()->GetDof(), BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
      dt(-1.0)
{
    linear_solver.iterative_mode = false;
    linear_solver.SetRelTol(1e-9);
    linear_solver.SetAbsTol(0.0);
    linear_solver.SetMaxIter(100);
    linear_solver.SetPrintLevel(0);
    linear_solver.SetPreconditioner(prec);
}

void DG_Solver::SetTimeStep(real_t dt_)
{
    if (dt_ != dt)
    {
        dt = dt_;
        A = K;
        A *= -dt;
        A += M;
        linear_solver.SetOperator(A);
    }
}

void DG_Solver::SetOperator(const Operator &op)
{
    linear_solver.SetOperator(op);
}

void DG_Solver::Mult(const Vector &x, Vector &y) const
{
    linear_solver.Mult(x, y);
}

FE_Evolution::FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_)
    : TimeDependentOperator(M_.FESpace()->GetTrueVSize()), M(M_), K(K_), b(b_), z(height)
{
    Array<int> ess_tdof_list;
    if (M.GetAssemblyLevel() == AssemblyLevel::LEGACY)
    {
        M_prec = new DSmoother(M.SpMat());
        M_solver.SetOperator(M.SpMat());
        dg_solver = new DG_Solver(M.SpMat(), K.SpMat(), *M.FESpace());
    }
    else
    {
        M_prec = new OperatorJacobiSmoother(M, ess_tdof_list);
        M_solver.SetOperator(M);
        dg_solver = NULL;
    }

    M_solver.SetPreconditioner(*M_prec);
    M_solver.iterative_mode = false;
    M_solver.SetRelTol(1e-9);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
    K.Mult(x, z);
    z += b;
    M_solver.Mult(z, y);
}

void FE_Evolution::ImplicitSolve(const real_t dt, const Vector &x, Vector &k)
{
    MFEM_VERIFY(dg_solver != NULL, "Implicit time integration is not supported with partial assembly");
    K.Mult(x, z);
    z += b;
    dg_solver->SetTimeStep(dt);
    dg_solver->Mult(z, k);
}

FE_Evolution::~FE_Evolution()
{
    delete M_prec;
    delete dg_solver;
}
