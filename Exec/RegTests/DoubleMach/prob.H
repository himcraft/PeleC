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
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  if ((x > prob_parm.x1) && (y < (x-prob_parm.x1)*1.7320508075688772)) {
    state(i, j, k, URHO) = prob_parm.rho0;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rho0 * prob_parm.e0;
    state(i, j, k, UEDEN) = prob_parm.rho0 * prob_parm.e0;
    state(i, j, k, UTEMP) = prob_parm.T0;
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, UFS + n) = prob_parm.rho0 * prob_parm.massfrac[n];

  } else {
    state(i, j, k, URHO) = prob_parm.rho1;
    state(i, j, k, UMX) = prob_parm.vx * prob_parm.rho1;
    state(i, j, k, UMY) = prob_parm.vy * prob_parm.rho1;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rho1 * prob_parm.e1;
    state(i, j, k, UEDEN) = prob_parm.rho1 * (prob_parm.e1 + 
		    0.5 * (prob_parm.vx * prob_parm.vx + 
			    prob_parm.vy * prob_parm.vy));
    state(i, j, k, UTEMP) = prob_parm.T1;
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, UFS + n) = prob_parm.rho1 * prob_parm.massfrac[n];
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
  const amrex::Real time,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  amrex::Real headx = time * 11.547005383792516 + prob_parm.x1 + 0.5773502691896257; 
  amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, eden = 0.0, T = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  if ((idir == 0) && (sgn == 1)) { 
    T = prob_parm.T1;
    u = prob_parm.vx * prob_parm.rho1;
    v = prob_parm.vy * prob_parm.rho1;
    rho = prob_parm.rho1;
    eint = prob_parm.rho1 * prob_parm.e1;
    eden = prob_parm.rho1 * (prob_parm.e1 +
		    0.5 * (prob_parm.vx * prob_parm.vx +
			    prob_parm.vy * prob_parm.vy));
  }
  if ((idir == 1) && (sgn == -1)) {
    if (x[0] >= headx){
      T = prob_parm.T0;
      u = 0.0;
      v = 0.0;
      rho = prob_parm.rho0;
      eint = prob_parm.rho0 * prob_parm.e0;
      eden = prob_parm.rho0 * prob_parm.e0;
    } else {
      T = prob_parm.T1;
      u = prob_parm.vx * prob_parm.rho1;
      v = prob_parm.vy * prob_parm.rho1;
      rho = prob_parm.rho1;
      eint = prob_parm.rho1 * prob_parm.e1;
      eden = prob_parm.rho1 * (prob_parm.e1 +
  		     0.5 * (prob_parm.vx * prob_parm.vx +
			     prob_parm.vy * prob_parm.vy));
    }
  } 
  if((idir == 1) && (sgn == 1)) {
    if (x[0] < prob_parm.x1) {
      T = prob_parm.T1;
      u = prob_parm.vx * prob_parm.rho1;
      v = prob_parm.vy * prob_parm.rho1;
      rho = prob_parm.rho1;
      eint = prob_parm.rho1 * prob_parm.e1;
      eden = prob_parm.rho1 * (prob_parm.e1 +
		     0.5 * (prob_parm.vx * prob_parm.vx +
			     prob_parm.vy * prob_parm.vy));
    } else {
      T = s_int[UTEMP];
      u = s_int[UMX];
      v = -s_int[UMY];
      rho = s_int[URHO];
      eint = s_int[UEINT];
      eden = s_int[UEDEN];
    }
  }
  s_ext[URHO] = rho;
  s_ext[UMX] = u;
  s_ext[UMY] = v;
  s_ext[UEINT] = eint;
  s_ext[UEDEN] = eden;
  s_ext[UTEMP] = T;  
  for (int n=0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * prob_parm.massfrac[n];}
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
