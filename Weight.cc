#include "Weight.h"


double calc_weight(const bsim::Decay& decay,
              const double det_angle,
              const double energy_ratio,
              const double parent_energy,
              const double nu_energy,
              const double gamma,
              const ROOT::RVec<double>& rr
              )
{
    double wght = (det_angle * decay.nimpwt * (energy_ratio * energy_ratio));

    if (std::abs(decay.ptype) != 13) {
        return wght;
    }

    //boost new neutrino to mu decay cm
    double beta[3];
    double p_nu[3]; //nu momentum
    beta[0] = decay.pdpx / parent_energy;
    beta[1] = decay.pdpy / parent_energy;
    beta[2] = decay.pdpz / parent_energy;

    const double rr_mag = std::sqrt(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]);

    p_nu[0] = rr[0] * nu_energy / rr_mag;
    p_nu[1] = rr[1] * nu_energy / rr_mag;
    p_nu[2] = rr[2] * nu_energy / rr_mag;

    double partial = gamma * (beta[0] * p_nu[0] + beta[1] * p_nu[1] + beta[2] * p_nu[2]);
    partial = nu_energy - partial / (gamma + 1.);
    double p_dcm_nu[4];
    for (int i = 0; i < 3; i++)
      p_dcm_nu[i] = p_nu[i] - beta[i] * gamma * partial;
    p_dcm_nu[3] = 0.;
    for (int i = 0; i < 3; i++)
      p_dcm_nu[3] += p_dcm_nu[i] * p_dcm_nu[i];
    p_dcm_nu[3] = std::sqrt(p_dcm_nu[3]);

    //boost parent of mu to mu production cm
    // gamma = m_ppenergy / m_parent_mass;
    beta[0] = decay.ppdxdz * decay.pppz / decay.ppenergy;
    beta[1] = decay.ppdydz * decay.pppz / decay.ppenergy;
    beta[2] = decay.pppz / decay.ppenergy;
    partial = gamma * (beta[0] * decay.muparpx + beta[1] * decay.muparpy +
                       beta[2] * decay.muparpz);
    partial = decay.mupare - partial / (gamma + 1.);
    double p_pcm_mp[4];
    p_pcm_mp[0] = decay.muparpx - beta[0] * gamma * partial;
    p_pcm_mp[1] = decay.muparpy - beta[1] * gamma * partial;
    p_pcm_mp[2] = decay.muparpz - beta[2] * gamma * partial;
    p_pcm_mp[3] = 0.;
    for (int i = 0; i < 3; i++)
      p_pcm_mp[3] += p_pcm_mp[i] * p_pcm_mp[i];
    p_pcm_mp[3] = std::sqrt(p_pcm_mp[3]);

    double wt_ratio = 1.;
    //have to check p_pcm_mp
    //it can be 0 if mupar..=0. (I guess muons created in target??)
    if (p_pcm_mp[3] != 0.) {
      //calc new decay angle w.r.t. (anti)spin direction
      double costh =
        (p_dcm_nu[0] * p_pcm_mp[0] + p_dcm_nu[1] * p_pcm_mp[1] + p_dcm_nu[2] * p_pcm_mp[2]) /
        (p_dcm_nu[3] * p_pcm_mp[3]);

      if (costh > 1.) costh = 1.;
      else if (costh < -1.) costh = -1.;

      //calc relative weight due to angle difference
      auto const nu_type = decay.ntype;
      if (std::abs(nu_type) == 12) { wt_ratio = 1. - costh; }
      else if (std::abs(nu_type) == 14) {
        constexpr double mumass = 0.13957039;
        double xnu = 2. * decay.necm / mumass;
        wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
      }
      else if (std::abs(nu_type) == 16) {
        constexpr double taumass = 1.77686;
        double xnu = 2. * decay.necm / taumass;
        wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
        std::cout << "calculating weight for tau neutrino; this may not be correct" << std::endl;
      }
      else {
        std::cout << "eventRates:: Bad neutrino type = " << nu_type << std::endl;
      }
    }
    return wght * wt_ratio;
}
