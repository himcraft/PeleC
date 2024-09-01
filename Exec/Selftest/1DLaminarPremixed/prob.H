#ifndef PROB_H
#define PROB_H

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
  //const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  amrex::Real u[3] = {0.0};
  if (x <= 0.1){
    state(i, j, k, URHO) = prob_parm.rho_i;
    state(i, j, k, UMX) = prob_parm.rho_i * u[0];
    state(i, j, k, UMY) = prob_parm.rho_i * u[1];
    state(i, j, k, UMZ) = prob_parm.rho_i * u[2];
    state(i, j, k, UEINT) = prob_parm.rho_i * prob_parm.e_i;
    state(i, j, k, UEDEN) = prob_parm.rho_i * (prob_parm.e_i +
		    0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    state(i, j, k, UTEMP) = prob_parm.T_i;
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, UFS + n) = prob_parm.rho_i * prob_parm.massfrac[n];

  } else {
    state(i, j, k, URHO) = prob_parm.rho;
    state(i, j, k, UMX) = prob_parm.rho_i * u[0];
    state(i, j, k, UMY) = prob_parm.rho_i * u[1];
    state(i, j, k, UMZ) = prob_parm.rho_i * u[2];
    state(i, j, k, UEINT) = prob_parm.rho * prob_parm.e;
    state(i, j, k, UEDEN) = prob_parm.rho_i * (prob_parm.e_i +
		    0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    state(i, j, k, UTEMP) = prob_parm.T;
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, UFS + n) = prob_parm.rho * prob_parm.massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, eden = 0.0, T = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == 1) { 
    T = prob_parm.T_i;
    u = 0.0;
    v = 0.0;
    w = 0.0;
    rho = prob_parm.rho_i;
    eint = prob_parm.rho_i * prob_parm.e_i;
    eden = prob_parm.rho_i * (prob_parm.e_i +
		    0.5 * (u * u + v * v + w * w));
  } else {
    T = prob_parm.T;
    u = 0.0;
    v = 0.0;
    w = 0.0;
    rho = prob_parm.rho;
    eint = prob_parm.rho * prob_parm.e;
    eden = prob_parm.rho_i * (prob_parm.e_i +
		    0.5 * (u * u + v * v + w * w));
  }
  s_ext[URHO] = rho;
  s_ext[UMX] = u;
  s_ext[UMY] = v;
  s_ext[UMZ] = w;
  s_ext[UEINT] = eint;
  s_ext[UEDEN] = eden;
  s_ext[UTEMP] = T;  
  for (int n=0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * prob_parm.massfrac[n];}
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
