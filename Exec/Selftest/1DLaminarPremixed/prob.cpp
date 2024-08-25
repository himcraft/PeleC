#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p0", PeleC::h_prob_parm_device->p);
    pp.query("rho0", PeleC::h_prob_parm_device->rho);
    pp.query("rho1", PeleC::h_prob_parm_device->rho_i);
    pp.query("T", PeleC::h_prob_parm_device->T);
    pp.query("T_i", PeleC::h_prob_parm_device->T_i);
    pp.query("e", PeleC::h_prob_parm_device->e);
    pp.query("e_i", PeleC::h_prob_parm_device->e_i);
    pp.query("Y_init_O2", PeleC::h_prob_parm_device->Y_init_O2);
    pp.query("Y_init_H2", PeleC::h_prob_parm_device->Y_init_H2);
    pp.query("Y_init_N2", PeleC::h_prob_parm_device->Y_init_N2);
  }

  PeleC::h_prob_parm_device->massfrac[H2_ID] = 
	  PeleC::h_prob_parm_device->Y_init_H2;
  PeleC::h_prob_parm_device->massfrac[O2_ID] =
	  PeleC::h_prob_parm_device->Y_init_O2;
  PeleC::h_prob_parm_device->massfrac[N2_ID] =
	  PeleC::h_prob_parm_device->Y_init_N2;

  auto eos = pele::physics::PhysicsType::eos();
  //amrex::Real cp = 0.0;
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->T,
    PeleC::h_prob_parm_device->rho,PeleC::h_prob_parm_device->e);

  eos.PYT2RE(
    PeleC::h_prob_parm_device->p,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->T_i,
    PeleC::h_prob_parm_device->rho_i,PeleC::h_prob_parm_device->e_i);

  /*eos.TY2Cp(
    PeleC::h_prob_parm_device->T, PeleC::h_prob_parm_device->massfrac.begin(),
    cp);

  auto& trans_parm = PeleC::trans_parms.host_parm();
  amrex::Real Pr = 0.7;
  trans_parm.const_conductivity = trans_parm.const_viscosity * cp / Pr;
  PeleC::trans_parms.sync_to_device();*/
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
