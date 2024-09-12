#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"

#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  amrex::Real u[3] = {0.0};

  // Set the states
  if (x <= 0.2) {
    state(i, j, k, URHO) = prob_parm.rho_l;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEDEN) = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l *
                                                 prob_parm.u_l *
                                                 prob_parm.u_l;
    state(i, j, k, UEINT) = prob_parm.rhoe_l;
    state(i, j, k, UTEMP) = prob_parm.T_l;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_l * prob_parm.massfrac_l[n];
    }
  } else {
    state(i, j, k, URHO) = prob_parm.rho_r;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEDEN) = prob_parm.rhoe_r + 0.5 * prob_parm.rho_r *
                                                 prob_parm.u_r *
                                                 prob_parm.u_r;
    state(i, j, k, UEINT) = prob_parm.rhoe_r;
    state(i, j, k, UTEMP) = prob_parm.T_r;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_r * prob_parm.massfrac_r[n];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.massfrac_r[n];
  }
  auto eos = pele::physics::PhysicsType::eos();

  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
  const amrex::Real rho_int = s_int[URHO];
  const amrex::Real u_int = s_int[UMX] / rho_int;
  const amrex::Real v_int = s_int[UMY] / rho_int;
  const amrex::Real w_int = s_int[UMZ] / rho_int;
  const amrex::Real T_int = s_int[UTEMP];
  amrex::Real p_int = 0.0, cs_int = 0.0;
  eos.RTY2P(rho_int, T_int, massfrac, p_int);
  eos.RTY2Cs(rho_int, T_int, massfrac, cs_int);

  const amrex::Real xb = prob_hi[0];
  const amrex::Real pb = prob_parm.p_r;
  const amrex::Real rhob =
    s_int[URHO] + (pb - p_int) / (cs_int * cs_int);
  const amrex::Real ub = u_int + (p_int - pb) / (rho_int * cs_int);
  const amrex::Real vb = v_int;
  const amrex::Real wb = w_int;

  rho = (rhob - rho_int) / (xb - xd) * (x[0] -xd) + rho_int;
  const amrex::Real p = (pb - p_int) / (xb - xd) * (x[0] - xd) + p_int;

  eos.RYP2E(rho, massfrac, p, eint);
  eos.EY2T(eint, massfrac, T);

  u = (ub - u_int) / (xb - xd) * (x[0] - xd) + u_int;
  v = (vb - v_int) / (xb - xd) * (x[0] - xd) + v_int;
  w = (wb - w_int) / (xb - xd) * (x[0] - xd) + w_int;

  s_ext[URHO] = rho;
  s_ext[UMX] = rho * u;
  s_ext[UMY] = rho * v;
  s_ext[UMZ] = rho * w;
  s_ext[UEINT] = rho * eint;
  s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
  s_ext[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * prob_parm.massfrac_r[n];
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
