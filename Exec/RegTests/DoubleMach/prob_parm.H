#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.0;
  amrex::Real p1 = 116.5;
  amrex::Real rho0 = 1.4;
  amrex::Real rho1 = 8.0;
  amrex::Real T0 = 300.0;
  amrex::Real T1 = 300.0;
  amrex::Real e0 = 1.0;
  amrex::Real e1 = 1.0;
  amrex::Real x1 = 0.2;
  amrex::Real vx = 7.144709581221618;
  amrex::Real vy = -4.125;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
