#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1e5;
  amrex::Real p1 = 1.165e7;
  amrex::Real rho0 = 0.014;
  amrex::Real rho1 = 0.08;
  amrex::Real T0 = 300.0;
  amrex::Real T1 = 300.0;
  amrex::Real e0 = 1.0;
  amrex::Real e1 = 1.0;
  amrex::Real x1 = 0.2;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
