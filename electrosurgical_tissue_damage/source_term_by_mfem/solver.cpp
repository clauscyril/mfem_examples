#include "solver.hpp"
#include <iostream>
#include <cmath>


using namespace mfem;


void compute_Q(const Mesh &mesh, FiniteElementSpace *fespace, GridFunction &q){
    int affichage = 1;
    std::cout << "Computing q using fem" << std::endl;

    // On précise quels noeuds sont des noeuds de Dirichlet
    Array<int> ess_tdof_list;
    for (int i = 0; i < mesh.GetNV(); i++) {
        const double *u = mesh.GetVertex(i);
        if ((u[0] + 0.35)*(u[0]+0.35) + (u[1])*(u[1]) <= 0.05*0.05 or (u[0] - 0.35)*(u[0]-0.35) + (u[1])*(u[1]) <= 0.05*0.05) {  // Vérifie si y = 0
            ess_tdof_list.Append(i);
        } 
    }

    GridFunction v(fespace);
    v = 0.f;

    for (int i = 0; i < ess_tdof_list.Size(); i++){
        const double *u = mesh.GetVertex(ess_tdof_list[i]);
        if (u[0] > 0) {
            v(ess_tdof_list[i]) = 1.f;
        } else {
            v(ess_tdof_list[i]) = -1.f;
        }
    }


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


    GradNormSquaredCoefficient grad_v_squared(v);

    q.ProjectCoefficient(grad_v_squared);

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
            sol_sock_q << "solution\n" << mesh << q
                        << "window_title 'Solution: q'" 
                        << "pause\n" << "keys c\n" << std::flush;
    }
        
}


GradNormSquaredCoefficient::GradNormSquaredCoefficient(const GridFunction &gf_)
    : gf(gf_), grad_coeff(&gf_),
      grad(gf_.FESpace()->GetMesh()->SpaceDimension()) {}

real_t GradNormSquaredCoefficient::Eval(ElementTransformation &T,
                                        const IntegrationPoint &ip)
{
    grad_coeff.Eval(grad, T, ip);
    return pow(grad.Norml2(),2);
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
