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
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p_l", PeleC::h_prob_parm_device->p_l);
    pp.query("T_l", PeleC::h_prob_parm_device->T_l);
    pp.query("p_r", PeleC::h_prob_parm_device->p_r);
    pp.query("T_r", PeleC::h_prob_parm_device->T_r);
    pp.query("EquivRatio", PeleC::h_prob_parm_device->ER);
  }

  PeleC::h_prob_parm_device->massfrac_r[H2_ID] = PeleC::h_prob_parm_device->ER /
	  (34.3196 + PeleC::h_prob_parm_device->ER);
  PeleC::h_prob_parm_device->massfrac_r[O2_ID] = 7.9964668 / 
	  (34.3196 + PeleC::h_prob_parm_device->ER); 
  PeleC::h_prob_parm_device->massfrac_r[N2_ID] = 1.0 - 
	  PeleC::h_prob_parm_device->massfrac_r[H2_ID] - PeleC::h_prob_parm_device->massfrac_r[O2_ID];
  PeleC::h_prob_parm_device->massfrac_l[N2_ID] = 1.0;

  amrex::Real e_l;
  amrex::Real e_r /*, cs, cp */;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_l, PeleC::h_prob_parm_device->massfrac_l.begin(),
    PeleC::h_prob_parm_device->T_l, PeleC::h_prob_parm_device->rho_l, e_l);
  PeleC::h_prob_parm_device->rhoe_l = PeleC::h_prob_parm_device->rho_l * e_l;
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_r, PeleC::h_prob_parm_device->massfrac_r.begin(),
    PeleC::h_prob_parm_device->T_r, PeleC::h_prob_parm_device->rho_r, e_r);
  PeleC::h_prob_parm_device->rhoe_r = PeleC::h_prob_parm_device->rho_r * e_r;
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
