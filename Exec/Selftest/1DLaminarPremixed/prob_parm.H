#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p = 1013250.0;
  amrex::Real rho = 1.4;
  amrex::Real rho_i = 8.0;
  amrex::Real T = 300.0;
  amrex::Real T_i = 300.0;
  amrex::Real e = 1.0;
  amrex::Real e_i = 1.0;
  amrex::Real Y_init_H2 = 0.06;
  amrex::Real Y_init_O2 = 0.5;
  amrex::Real Y_init_N2 = 0.44;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
