#include "solver.hpp"

using namespace mfem;

real_t q_func(const Vector &x, real_t &t)
{
    if (t < 1000)
        if ((x(0) < -0.45) || (x(0) > 0.45))
            return 1.;
        else 
            return 0.;
    else 
        return 0.;
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
