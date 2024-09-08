#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_l = 1.0;   // left pressure (erg/cc)
  amrex::Real u_l = 0.0;   // left velocity (cm/s)
  amrex::Real rho_l = 0.0; // left density (g/cc)
  amrex::Real rhoe_l;
  amrex::Real T_l = 1.0;
  amrex::Real p_r = 0.1;     // right pressure (erg/cc)
  amrex::Real u_r = 0.0;     // right velocity (cm/s)
  amrex::Real rho_r = 0.125; // right density (g/cc)
  amrex::Real rhoe_r;
  amrex::Real T_r = 1.0;
  amrex::Real ER = 1.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_r = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_l = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
