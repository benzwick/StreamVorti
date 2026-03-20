/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * Space-time integrators implementation.
 *
 * Implements SUPG/PSPG/LSIC stabilized finite element integrators
 * for the incompressible Navier-Stokes equations in space-time.
 */

#include <StreamVorti/spacetime/st_integrators.hpp>
#include <cmath>
#include <algorithm>

namespace StreamVorti {
namespace SpaceTime {

// =========================================================================
// Stabilization parameter computations
// =========================================================================

double ComputeTauSUGN(const mfem::Vector &velocity,
                      double h_elem, double dt, double nu,
                      const StabilizationParams &params)
{
    double u_norm = velocity.Norml2();

    // Three time scales (Tezduyar 2003):
    // τ_1: advective time scale
    double inv_tau1 = 0.0;
    if (u_norm > 1e-14)
    {
        inv_tau1 = params.c_1 * u_norm / h_elem;
    }

    // τ_2: temporal time scale
    double inv_tau2 = params.c_2 / dt;

    // τ_3: diffusive time scale
    double inv_tau3 = 0.0;
    if (nu > 1e-14)
    {
        inv_tau3 = params.c_3 * nu / (h_elem * h_elem);
    }

    // Combined: 1/τ = sqrt(1/τ₁² + 1/τ₂² + 1/τ₃²)
    double inv_tau_sq = inv_tau1 * inv_tau1 +
                        inv_tau2 * inv_tau2 +
                        inv_tau3 * inv_tau3;

    return (inv_tau_sq > 1e-30) ? 1.0 / std::sqrt(inv_tau_sq) : 0.0;
}

double ComputeTauLSIC(const mfem::Vector &velocity,
                      double h_elem, double tau_sugn)
{
    double u_norm = velocity.Norml2();

    // τ_LSIC = h |u| / 2 (simplified)
    // Or equivalently: τ_LSIC = h² / τ_SUGN (Behr's form)
    if (tau_sugn > 1e-30)
    {
        return h_elem * h_elem / tau_sugn;
    }
    return h_elem * u_norm * 0.5;
}

double ComputeElementLength(const mfem::Vector &velocity,
                            mfem::ElementTransformation &Trans)
{
    // Element length in the flow direction using the inverse Jacobian metric
    // h = 2 / ||J^{-T} s||  where s = u/|u| is the unit flow direction

    double u_norm = velocity.Norml2();
    if (u_norm < 1e-14)
    {
        // No flow direction: use isotropic element size
        return std::pow(Trans.Weight(), 1.0 / Trans.GetDimension());
    }

    const mfem::IntegrationPoint &ip = Trans.GetIntPoint();
    const mfem::DenseMatrix &Jinv = Trans.InverseJacobian();

    int dim = velocity.Size();
    mfem::Vector s(dim);
    s = velocity;
    s /= u_norm;

    // r = J^{-T} s
    mfem::Vector r(Jinv.Height());
    Jinv.MultTranspose(s, r);

    double r_norm = r.Norml2();
    return (r_norm > 1e-14) ? 2.0 / r_norm : Trans.Weight();
}

// =========================================================================
// STConvectionIntegrator
// =========================================================================

void STConvectionIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el,
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd = el.GetDof();
    int dim = el.GetDim(); // space-time dimension
    int sdim = dim - 1;    // spatial dimension

    elmat.SetSize(nd * sdim, nd * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 1);

    mfem::Vector shape(nd);
    mfem::DenseMatrix dshape(nd, dim);
    mfem::DenseMatrix dshape_phys(nd, dim);
    mfem::Vector vel(dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        el.CalcShape(ip, shape);
        el.CalcDShape(ip, dshape);

        // Transform derivatives to physical space
        mfem::Mult(dshape, Trans.InverseJacobian(), dshape_phys);

        // Get velocity at this point (spatial components + temporal = 1)
        velocity_coeff_.Eval(vel, Trans, ip);

        // Space-time advection velocity û = (u₁, u₂, ..., 1)
        // The last component (time direction) is 1
        mfem::Vector st_vel(dim);
        for (int d = 0; d < sdim; ++d)
        {
            st_vel(d) = vel(d);
        }
        st_vel(dim - 1) = 1.0; // time direction velocity = 1

        // Compute û · ∇N for each shape function
        mfem::Vector u_dot_grad_N(nd);
        dshape_phys.Mult(st_vel, u_dot_grad_N);

        // Assemble: ρ (w_i, û·∇N_j u_j) for each velocity component
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd; ++i)
            {
                for (int j = 0; j < nd; ++j)
                {
                    elmat(d * nd + i, d * nd + j) +=
                        rho_ * w * shape(i) * u_dot_grad_N(j);
                }
            }
        }
    }
}

// =========================================================================
// STDiffusionIntegrator
// =========================================================================

void STDiffusionIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el,
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd = el.GetDof();
    int dim = el.GetDim();
    int sdim = dim - 1; // spatial dimension

    elmat.SetSize(nd * sdim, nd * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());

    mfem::DenseMatrix dshape(nd, dim);
    mfem::DenseMatrix dshape_phys(nd, dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        el.CalcDShape(ip, dshape);
        mfem::Mult(dshape, Trans.InverseJacobian(), dshape_phys);

        // Only spatial derivatives contribute to viscous term
        // (∇_x w, μ ∇_x u) - no time derivatives
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd; ++i)
            {
                for (int j = 0; j < nd; ++j)
                {
                    double val = 0.0;
                    for (int k = 0; k < sdim; ++k) // spatial dims only
                    {
                        val += dshape_phys(i, k) * dshape_phys(j, k);
                    }
                    elmat(d * nd + i, d * nd + j) += mu_ * w * val;
                }
            }
        }
    }
}

// =========================================================================
// STPressureGradientIntegrator
// =========================================================================

void STPressureGradientIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &trial_el,  // pressure
    const mfem::FiniteElement &test_el,   // velocity
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd_v = test_el.GetDof();
    int nd_p = trial_el.GetDof();
    int dim = test_el.GetDim();
    int sdim = dim - 1;

    elmat.SetSize(nd_v * sdim, nd_p);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(test_el.GetGeomType(),
                           test_el.GetOrder() + trial_el.GetOrder());

    mfem::Vector shape_p(nd_p);
    mfem::DenseMatrix dshape_v(nd_v, dim);
    mfem::DenseMatrix dshape_v_phys(nd_v, dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        trial_el.CalcShape(ip, shape_p);
        test_el.CalcDShape(ip, dshape_v);
        mfem::Mult(dshape_v, Trans.InverseJacobian(), dshape_v_phys);

        // -(∇_x · w, p) = -Σ_d (∂w_d/∂x_d, p)
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd_v; ++i)
            {
                for (int j = 0; j < nd_p; ++j)
                {
                    elmat(d * nd_v + i, j) -=
                        w * dshape_v_phys(i, d) * shape_p(j);
                }
            }
        }
    }
}

// =========================================================================
// STDivergenceIntegrator
// =========================================================================

void STDivergenceIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &trial_el,  // velocity
    const mfem::FiniteElement &test_el,   // pressure
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd_v = trial_el.GetDof();
    int nd_p = test_el.GetDof();
    int dim = trial_el.GetDim();
    int sdim = dim - 1;

    elmat.SetSize(nd_p, nd_v * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(trial_el.GetGeomType(),
                           trial_el.GetOrder() + test_el.GetOrder());

    mfem::Vector shape_p(nd_p);
    mfem::DenseMatrix dshape_v(nd_v, dim);
    mfem::DenseMatrix dshape_v_phys(nd_v, dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        test_el.CalcShape(ip, shape_p);
        trial_el.CalcDShape(ip, dshape_v);
        mfem::Mult(dshape_v, Trans.InverseJacobian(), dshape_v_phys);

        // (q, ∇_x · u) = Σ_d (q, ∂u_d/∂x_d)
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd_p; ++i)
            {
                for (int j = 0; j < nd_v; ++j)
                {
                    elmat(i, d * nd_v + j) +=
                        w * shape_p(i) * dshape_v_phys(j, d);
                }
            }
        }
    }
}

// =========================================================================
// STSUPGIntegrator
// =========================================================================

void STSUPGIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el,
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd = el.GetDof();
    int dim = el.GetDim();
    int sdim = dim - 1;

    elmat.SetSize(nd * sdim, nd * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 1);

    mfem::Vector shape(nd);
    mfem::DenseMatrix dshape(nd, dim);
    mfem::DenseMatrix dshape_phys(nd, dim);
    mfem::Vector vel(dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        el.CalcShape(ip, shape);
        el.CalcDShape(ip, dshape);
        mfem::Mult(dshape, Trans.InverseJacobian(), dshape_phys);

        velocity_coeff_.Eval(vel, Trans, ip);

        // Space-time advection velocity
        mfem::Vector st_vel(dim);
        for (int d = 0; d < sdim; ++d)
        {
            st_vel(d) = vel(d);
        }
        st_vel(dim - 1) = 1.0;

        // Compute element length and stabilization parameter
        double h_e = ComputeElementLength(vel, Trans);
        double dt_local = std::pow(Trans.Weight(), 1.0 / dim);
        double tau = ComputeTauSUGN(vel, h_e, dt_local, mu_ / rho_, params_);

        // Compute û · ∇N for each shape function
        mfem::Vector u_dot_grad_N(nd);
        dshape_phys.Mult(st_vel, u_dot_grad_N);

        // SUPG: τ (ρ û·∇w) · (ρ û·∇u)
        // This is the advection-advection part of SUPG
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd; ++i)
            {
                for (int j = 0; j < nd; ++j)
                {
                    elmat(d * nd + i, d * nd + j) +=
                        tau * rho_ * rho_ * w *
                        u_dot_grad_N(i) * u_dot_grad_N(j);
                }
            }
        }

        // SUPG: τ (ρ û·∇w) · (-μ ∇²u)
        // This is the advection-diffusion cross term
        // Simplified: we include only the main diagonal contribution
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd; ++i)
            {
                for (int j = 0; j < nd; ++j)
                {
                    double laplacian_j = 0.0;
                    // Approximate ∇²N_j using second derivatives
                    // (not available directly; this term is often dropped
                    //  for linear elements as ∇²N = 0)
                    // For higher-order elements, this could use
                    // CalcHessian if available.

                    elmat(d * nd + i, d * nd + j) -=
                        tau * rho_ * mu_ * w *
                        u_dot_grad_N(i) * laplacian_j;
                }
            }
        }
    }
}

// =========================================================================
// STPSPGIntegrator
// =========================================================================

void STPSPGIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &trial_el,  // velocity
    const mfem::FiniteElement &test_el,   // pressure
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd_v = trial_el.GetDof();
    int nd_p = test_el.GetDof();
    int dim = trial_el.GetDim();
    int sdim = dim - 1;

    elmat.SetSize(nd_p, nd_v * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(trial_el.GetGeomType(),
                           trial_el.GetOrder() + test_el.GetOrder() + 1);

    mfem::Vector shape_v(nd_v);
    mfem::DenseMatrix dshape_v(nd_v, dim);
    mfem::DenseMatrix dshape_v_phys(nd_v, dim);
    mfem::DenseMatrix dshape_p(nd_p, dim);
    mfem::DenseMatrix dshape_p_phys(nd_p, dim);
    mfem::Vector vel(dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        trial_el.CalcShape(ip, shape_v);
        trial_el.CalcDShape(ip, dshape_v);
        mfem::Mult(dshape_v, Trans.InverseJacobian(), dshape_v_phys);

        test_el.CalcDShape(ip, dshape_p);
        mfem::Mult(dshape_p, Trans.InverseJacobian(), dshape_p_phys);

        velocity_coeff_.Eval(vel, Trans, ip);

        mfem::Vector st_vel(dim);
        for (int d = 0; d < sdim; ++d)
        {
            st_vel(d) = vel(d);
        }
        st_vel(dim - 1) = 1.0;

        double h_e = ComputeElementLength(vel, Trans);
        double dt_local = std::pow(Trans.Weight(), 1.0 / dim);
        double tau = ComputeTauSUGN(vel, h_e, dt_local, mu_ / rho_, params_);

        // Compute û · ∇N_v for velocity shape functions
        mfem::Vector u_dot_grad_Nv(nd_v);
        dshape_v_phys.Mult(st_vel, u_dot_grad_Nv);

        // PSPG: τ/ρ (∇_x q) · (ρ û·∇u)
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd_p; ++i)
            {
                for (int j = 0; j < nd_v; ++j)
                {
                    elmat(i, d * nd_v + j) +=
                        tau * w * dshape_p_phys(i, d) *
                        rho_ * u_dot_grad_Nv(j);
                }
            }
        }
    }
}

// =========================================================================
// STLSICIntegrator
// =========================================================================

void STLSICIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el,
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd = el.GetDof();
    int dim = el.GetDim();
    int sdim = dim - 1;

    elmat.SetSize(nd * sdim, nd * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());

    mfem::DenseMatrix dshape(nd, dim);
    mfem::DenseMatrix dshape_phys(nd, dim);
    mfem::Vector vel(dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        el.CalcDShape(ip, dshape);
        mfem::Mult(dshape, Trans.InverseJacobian(), dshape_phys);

        velocity_coeff_.Eval(vel, Trans, ip);

        double h_e = ComputeElementLength(vel, Trans);
        double dt_local = std::pow(Trans.Weight(), 1.0 / dim);
        double tau_sugn = ComputeTauSUGN(vel, h_e, dt_local,
                                          mu_ / rho_, params_);
        double tau_lsic = ComputeTauLSIC(vel, h_e, tau_sugn);

        // LSIC: τ_LSIC (∇_x · w)(∇_x · u)
        for (int i = 0; i < nd; ++i)
        {
            for (int j = 0; j < nd; ++j)
            {
                // div(w_i) · div(u_j) across all velocity components
                for (int d1 = 0; d1 < sdim; ++d1)
                {
                    for (int d2 = 0; d2 < sdim; ++d2)
                    {
                        elmat(d1 * nd + i, d2 * nd + j) +=
                            tau_lsic * rho_ * w *
                            dshape_phys(i, d1) * dshape_phys(j, d2);
                    }
                }
            }
        }
    }
}

// =========================================================================
// STTemporalIntegrator
// =========================================================================

void STTemporalIntegrator::AssembleElementMatrix(
    const mfem::FiniteElement &el,
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat)
{
    int nd = el.GetDof();
    int dim = el.GetDim();
    int sdim = dim - 1;

    elmat.SetSize(nd * sdim, nd * sdim);
    elmat = 0.0;

    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());

    mfem::Vector shape(nd);
    mfem::DenseMatrix dshape(nd, dim);
    mfem::DenseMatrix dshape_phys(nd, dim);

    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
        const mfem::IntegrationPoint &ip = ir.IntPoint(q);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        el.CalcShape(ip, shape);
        el.CalcDShape(ip, dshape);
        mfem::Mult(dshape, Trans.InverseJacobian(), dshape_phys);

        // (w, ρ ∂u/∂t) where ∂/∂t is the last coordinate direction
        for (int d = 0; d < sdim; ++d)
        {
            for (int i = 0; i < nd; ++i)
            {
                for (int j = 0; j < nd; ++j)
                {
                    // N_i * ρ * ∂N_j/∂t
                    elmat(d * nd + i, d * nd + j) +=
                        rho_ * w * shape(i) * dshape_phys(j, dim - 1);
                }
            }
        }
    }
}

// =========================================================================
// STNavierStokesOperator
// =========================================================================

STNavierStokesOperator::STNavierStokesOperator(
    mfem::FiniteElementSpace &vel_fes,
    mfem::FiniteElementSpace &pres_fes,
    const SolverConfig &config)
    : mfem::Operator(vel_fes.GetTrueVSize() + pres_fes.GetTrueVSize()),
      vel_fes_(vel_fes),
      pres_fes_(pres_fes),
      config_(config),
      jacobian_(nullptr),
      J_uu_(nullptr),
      J_up_(nullptr),
      J_pu_(nullptr),
      J_pp_(nullptr),
      body_force_(nullptr),
      has_prev_solution_(false)
{
    // Set up block structure
    block_offsets_.SetSize(3);
    block_offsets_[0] = 0;
    block_offsets_[1] = vel_fes.GetTrueVSize();
    block_offsets_[2] = vel_fes.GetTrueVSize() + pres_fes.GetTrueVSize();
}

STNavierStokesOperator::~STNavierStokesOperator()
{
    delete jacobian_;
    delete J_uu_;
    delete J_up_;
    delete J_pu_;
    delete J_pp_;
}

void STNavierStokesOperator::SetState(const mfem::BlockVector &state)
{
    // Store current state for linearization
    // This is called before each Newton iteration
}

void STNavierStokesOperator::SetPreviousSolution(const mfem::Vector &u_prev)
{
    u_prev_ = u_prev;
    has_prev_solution_ = true;
}

void STNavierStokesOperator::Mult(const mfem::Vector &x,
                                   mfem::Vector &y) const
{
    mfem::BlockVector bx(const_cast<mfem::Vector &>(x).GetData(),
                         block_offsets_);
    mfem::BlockVector by(y.GetData(), block_offsets_);
    AssembleResidual(bx, by);
}

mfem::Operator &STNavierStokesOperator::GetGradient(
    const mfem::Vector &x) const
{
    mfem::BlockVector bx(const_cast<mfem::Vector &>(x).GetData(),
                         block_offsets_);
    AssembleJacobian(bx);
    return *jacobian_;
}

void STNavierStokesOperator::AssembleResidual(
    const mfem::BlockVector &state,
    mfem::BlockVector &residual) const
{
    residual = 0.0;

    const mfem::Vector &u = state.GetBlock(0);
    const mfem::Vector &p = state.GetBlock(1);

    int sdim = config_.spatial_dim;
    double rho = config_.density;
    double mu = config_.viscosity;

    // Create a velocity grid function for the current state
    mfem::GridFunction u_gf(&vel_fes_);
    u_gf = u;

    mfem::GridFunction p_gf(&pres_fes_);
    p_gf = p;

    // Compute residual element-by-element
    mfem::Vector &R_M = residual.GetBlock(0); // Momentum residual
    mfem::Vector &R_C = residual.GetBlock(1); // Continuity residual

    // Assemble using MFEM's BilinearForm/LinearForm infrastructure
    // For Newton's method, the residual R(U) = K(U)*U - F
    // where K is the assembled operator matrix

    // This is a placeholder for the full nonlinear assembly.
    // In practice, each integrator contributes to the residual.
    // The actual implementation uses element-level assembly
    // with the custom integrators defined above.

    if (jacobian_)
    {
        jacobian_->Mult(state, residual);
    }
}

void STNavierStokesOperator::AssembleJacobian(
    const mfem::BlockVector &state) const
{
    int sdim = config_.spatial_dim;
    double rho = config_.density;
    double mu = config_.viscosity;

    // Clean up previous Jacobian
    delete jacobian_;
    delete J_uu_;
    delete J_up_;
    delete J_pu_;
    delete J_pp_;

    // Create velocity coefficient from current state
    mfem::GridFunction u_gf(&vel_fes_);
    u_gf = state.GetBlock(0);

    mfem::VectorGridFunctionCoefficient vel_coeff(&u_gf);

    // ---- Assemble J_uu (momentum-velocity block) ----
    {
        mfem::BilinearForm a(&vel_fes_);

        // Temporal derivative
        a.AddDomainIntegrator(new STTemporalIntegrator(rho));

        // Convection (linearized: advection by current velocity)
        a.AddDomainIntegrator(new STConvectionIntegrator(vel_coeff, rho));

        // Diffusion
        a.AddDomainIntegrator(new STDiffusionIntegrator(mu));

        // SUPG stabilization
        if (config_.stabilization.enable_supg)
        {
            a.AddDomainIntegrator(
                new STSUPGIntegrator(vel_coeff, rho, mu,
                                     config_.stabilization));
        }

        // LSIC stabilization
        if (config_.stabilization.enable_lsic)
        {
            a.AddDomainIntegrator(
                new STLSICIntegrator(vel_coeff, rho, mu,
                                     config_.stabilization));
        }

        a.Assemble();
        a.Finalize();
        J_uu_ = new mfem::SparseMatrix(a.SpMat());
    }

    // ---- Assemble J_up (momentum-pressure block) ----
    {
        mfem::MixedBilinearForm b(&pres_fes_, &vel_fes_);
        b.AddDomainIntegrator(new STPressureGradientIntegrator());
        b.Assemble();
        b.Finalize();
        J_up_ = new mfem::SparseMatrix(b.SpMat());
    }

    // ---- Assemble J_pu (continuity-velocity block) ----
    {
        mfem::MixedBilinearForm c(&vel_fes_, &pres_fes_);
        c.AddDomainIntegrator(new STDivergenceIntegrator());

        // PSPG stabilization (pressure equation)
        if (config_.stabilization.enable_pspg)
        {
            c.AddDomainIntegrator(
                new STPSPGIntegrator(vel_coeff, rho, mu,
                                     config_.stabilization));
        }

        c.Assemble();
        c.Finalize();
        J_pu_ = new mfem::SparseMatrix(c.SpMat());
    }

    // ---- Assemble J_pp (pressure-pressure block from PSPG) ----
    {
        mfem::BilinearForm d(&pres_fes_);
        // PSPG contributes a pressure Laplacian-like term:
        // τ (∇q, ∇p)
        if (config_.stabilization.enable_pspg)
        {
            d.AddDomainIntegrator(new mfem::DiffusionIntegrator());
        }
        d.Assemble();
        d.Finalize();
        J_pp_ = new mfem::SparseMatrix(d.SpMat());
    }

    // Build block operator
    jacobian_ = new mfem::BlockOperator(block_offsets_);
    jacobian_->SetBlock(0, 0, J_uu_);
    jacobian_->SetBlock(0, 1, J_up_);
    jacobian_->SetBlock(1, 0, J_pu_);
    jacobian_->SetBlock(1, 1, J_pp_);
}

} // namespace SpaceTime
} // namespace StreamVorti
