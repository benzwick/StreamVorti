/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * Space-time finite element integrators for incompressible Navier-Stokes.
 *
 * Implements the SUPG/PSPG/LSIC stabilized formulation in (d+1)-dimensional
 * space-time domain following Tezduyar, Behr et al.
 *
 * The incompressible Navier-Stokes equations in space-time weak form:
 *
 *   Galerkin terms:
 *     (w, ρ ∂u/∂t) + (w, ρ u·∇u) + (ε(w), σ) + (q, ∇·u) = (w, f)
 *
 *   where σ = -pI + 2μ ε(u), ε(u) = (∇u + ∇u^T)/2
 *
 *   SUPG stabilization:
 *     Σ_e ∫_Ωe τ_SUPG (u·∇w + ∂w/∂t) · R_M dΩ
 *
 *   PSPG stabilization:
 *     Σ_e ∫_Ωe τ_PSPG ∇q · R_M dΩ
 *
 *   LSIC stabilization:
 *     Σ_e ∫_Ωe τ_LSIC (∇·w)(∇·u) dΩ
 *
 *   where R_M = ρ(∂u/∂t + u·∇u) + ∇p - μ∇²u - f  (momentum residual)
 *
 * References:
 *   [1] Tezduyar, Behr, Liou (1992). "A new strategy for finite element
 *       computations involving moving boundaries and interfaces."
 *   [2] Tezduyar (2003). "Stabilization parameters in SUPG and PSPG."
 *   [3] Bazilevs, Calo, Cottrell, Hughes, Reali, Scovazzi (2007).
 *       "Variational multiscale residual-based turbulence modeling."
 */

#ifndef STREAMVORTI_SPACETIME_INTEGRATORS_HPP
#define STREAMVORTI_SPACETIME_INTEGRATORS_HPP

#include "st_config.hpp"
#include "mfem.hpp"

namespace StreamVorti {
namespace SpaceTime {

/// Compute stabilization parameter τ_SUPG using Tezduyar's formulation.
///
/// τ_SUGN = 1 / sqrt( (2/Δt)² + (2|u|/h)² + (4ν/h²)² )
///
/// This combines three time scales:
///   τ_SUGN1 = h / (2|u|)       - convective time scale
///   τ_SUGN2 = Δt / 2           - temporal time scale
///   τ_SUGN3 = h² / (4ν)        - diffusive time scale
///
/// For space-time, h is the element length in the flow direction and
/// Δt is absorbed into the element metric.
double ComputeTauSUGN(const mfem::Vector &velocity,
                      double h_elem, double dt, double nu,
                      const StabilizationParams &params);

/// Compute LSIC stabilization parameter τ_LSIC
///
/// τ_LSIC = h * |u| / 2  (simplified form)
/// or τ_LSIC = h² / τ_SUGN  (Behr's form)
double ComputeTauLSIC(const mfem::Vector &velocity,
                      double h_elem, double tau_sugn);

/// Compute element length h in the direction of flow
/// Uses the metric tensor approach from Tezduyar (2003)
double ComputeElementLength(const mfem::Vector &velocity,
                            mfem::ElementTransformation &Trans);

// =========================================================================
// Block integrators for the velocity-pressure system
// =========================================================================

/// Space-time convection integrator for momentum equation.
///
/// Implements: (w, ρ u·∇u) in space-time, where ∂u/∂t is treated
/// as part of the convective derivative in the time direction.
///
/// In space-time coordinates x̂ = (x, t), the advection velocity
/// includes the temporal component: û = (u, 1).
class STConvectionIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STConvectionIntegrator(mfem::VectorCoefficient &velocity,
                           double rho = 1.0)
        : velocity_coeff_(velocity), rho_(rho) {}

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override;

private:
    mfem::VectorCoefficient &velocity_coeff_;
    double rho_;
};

/// Space-time diffusion integrator for viscous term.
///
/// Implements: (ε(w), 2μ ε(u)) = μ (∇w, ∇u + ∇u^T)
/// Only acts on spatial dimensions (not time).
class STDiffusionIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STDiffusionIntegrator(double mu) : mu_(mu) {}

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override;

private:
    double mu_;
};

/// Pressure gradient integrator: -(∇·w, p)
/// Acts only on spatial gradient (not time derivative).
class STPressureGradientIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STPressureGradientIntegrator() {}

    void AssembleElementMatrix2(const mfem::FiniteElement &trial_el,
                                const mfem::FiniteElement &test_el,
                                mfem::ElementTransformation &Trans,
                                mfem::DenseMatrix &elmat) override;
    // Single-element version (not used for mixed)
    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override
    { /* unused for mixed FE */ }
};

/// Divergence integrator: (q, ∇·u)
/// Acts only on spatial divergence.
class STDivergenceIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STDivergenceIntegrator() {}

    void AssembleElementMatrix2(const mfem::FiniteElement &trial_el,
                                const mfem::FiniteElement &test_el,
                                mfem::ElementTransformation &Trans,
                                mfem::DenseMatrix &elmat) override;

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override
    { /* unused for mixed FE */ }
};

/// SUPG stabilization integrator.
///
/// Adds: Σ_e ∫_Ωe τ_SUPG (ρ û·∇w) · (ρ û·∇u + ∇p - μ∇²u) dΩ
///
/// where û = (u, 1) is the space-time advection velocity.
class STSUPGIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STSUPGIntegrator(mfem::VectorCoefficient &velocity,
                     double rho, double mu,
                     const StabilizationParams &params)
        : velocity_coeff_(velocity), rho_(rho), mu_(mu), params_(params) {}

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override;

private:
    mfem::VectorCoefficient &velocity_coeff_;
    double rho_;
    double mu_;
    StabilizationParams params_;
};

/// PSPG stabilization integrator.
///
/// Adds: Σ_e ∫_Ωe τ_PSPG (∇q) · (ρ û·∇u + ∇p - μ∇²u) dΩ
///
/// This stabilizes the pressure equation for equal-order interpolation.
class STPSPGIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STPSPGIntegrator(mfem::VectorCoefficient &velocity,
                     double rho, double mu,
                     const StabilizationParams &params)
        : velocity_coeff_(velocity), rho_(rho), mu_(mu), params_(params) {}

    void AssembleElementMatrix2(const mfem::FiniteElement &trial_el,
                                const mfem::FiniteElement &test_el,
                                mfem::ElementTransformation &Trans,
                                mfem::DenseMatrix &elmat) override;

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override
    { /* unused for mixed FE */ }

private:
    mfem::VectorCoefficient &velocity_coeff_;
    double rho_;
    double mu_;
    StabilizationParams params_;
};

/// LSIC stabilization integrator.
///
/// Adds: Σ_e ∫_Ωe τ_LSIC (∇·w)(∇·u) dΩ
///
/// Provides additional stability for the incompressibility constraint.
class STLSICIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STLSICIntegrator(mfem::VectorCoefficient &velocity,
                     double rho, double mu,
                     const StabilizationParams &params)
        : velocity_coeff_(velocity), rho_(rho), mu_(mu), params_(params) {}

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override;

private:
    mfem::VectorCoefficient &velocity_coeff_;
    double rho_;
    double mu_;
    StabilizationParams params_;
};

/// Temporal derivative integrator for space-time: (w, ρ ∂u/∂t)
///
/// In the space-time setting, this is the mass matrix in the
/// time direction. For dG in time, this includes the jump term
/// at the bottom of the time slab.
class STTemporalIntegrator : public mfem::BilinearFormIntegrator
{
public:
    STTemporalIntegrator(double rho = 1.0) : rho_(rho) {}

    void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override;

private:
    double rho_;
};

// =========================================================================
// Nonlinear integrators for Newton linearization
// =========================================================================

/// Combined nonlinear space-time Navier-Stokes operator.
///
/// Assembles the full nonlinear residual:
///   R(U) = [R_M(u,p); R_C(u)] where U = (u, p)
///
///   R_M = ρ(∂u/∂t + u·∇u) + ∇p - μ∇²u - f   (momentum)
///   R_C = ∇·u                                   (continuity)
///
/// plus SUPG/PSPG/LSIC stabilization terms.
///
/// The Jacobian J = dR/dU is assembled for Newton's method.
class STNavierStokesOperator : public mfem::Operator
{
public:
    STNavierStokesOperator(mfem::FiniteElementSpace &vel_fes,
                           mfem::FiniteElementSpace &pres_fes,
                           const SolverConfig &config);

    ~STNavierStokesOperator();

    /// Set the current solution for linearization
    void SetState(const mfem::BlockVector &state);

    /// Set body force
    void SetBodyForce(mfem::VectorCoefficient *f) { body_force_ = f; }

    /// Set the solution from the previous time slab (for dG jump term)
    void SetPreviousSolution(const mfem::Vector &u_prev);

    /// Compute the residual R(U) at the current state
    void Mult(const mfem::Vector &x, mfem::Vector &y) const override;

    /// Get the Jacobian dR/dU
    mfem::Operator &GetGradient(const mfem::Vector &x) const override;

    /// Get block offsets for the velocity-pressure system
    const mfem::Array<int> &GetBlockOffsets() const { return block_offsets_; }

    /// Get the spatial dimension
    int GetSpatialDim() const { return config_.spatial_dim; }

private:
    void AssembleResidual(const mfem::BlockVector &state,
                          mfem::BlockVector &residual) const;

    void AssembleJacobian(const mfem::BlockVector &state) const;

    mfem::FiniteElementSpace &vel_fes_;
    mfem::FiniteElementSpace &pres_fes_;
    SolverConfig config_;

    mutable mfem::BlockOperator *jacobian_;
    mutable mfem::SparseMatrix *J_uu_;  ///< d(R_M)/du
    mutable mfem::SparseMatrix *J_up_;  ///< d(R_M)/dp
    mutable mfem::SparseMatrix *J_pu_;  ///< d(R_C)/du
    mutable mfem::SparseMatrix *J_pp_;  ///< d(R_C)/dp (from PSPG)

    mfem::Array<int> block_offsets_;
    mfem::VectorCoefficient *body_force_;
    mfem::Vector u_prev_;  ///< Solution from previous slab
    bool has_prev_solution_;
};

} // namespace SpaceTime
} // namespace StreamVorti

#endif // STREAMVORTI_SPACETIME_INTEGRATORS_HPP
