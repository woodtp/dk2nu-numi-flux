#include "Weight.h"

std::vector<double> calc_pT(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<double> pT(ancestors.size() - 1);

  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if (i == 0) {
      pT[i] = -9999.;
      continue;
    }
    else if (ancestors[i - 1].nucleus == 0 || ancestors[i].nucleus == 1000000000) {
      pT[i] = -9999.;
      continue;
    }
    auto const& prod = ancestors[i];

    auto const p_prod = std::sqrt(prod.startpx * prod.startpx + prod.startpy * prod.startpy +
                                  prod.startpz * prod.startpz);
    auto const p_inc = std::sqrt(prod.pprodpx * prod.pprodpx + prod.pprodpy * prod.pprodpy +
                                 prod.pprodpz * prod.pprodpz);

    double costh =
      (prod.startpx * prod.pprodpx + prod.startpy * prod.pprodpy + prod.startpz * prod.pprodpz) /
      (p_prod * p_inc);

    if (costh > 1.)
      costh = 1.;
    else if (costh < -1.)
      costh = -1.;

    const double sinth = std::sqrt(1. - costh * costh);

    pT[i] = p_prod * sinth;
  }

  return pT;
}

std::vector<double> calc_xF(const ROOT::RVec<bsim::Ancestor>& ancestors,
                            const ROOT::RVec<double>& ancestor_masses)
{
  std::vector<double> xF(ancestors.size() - 1);

  // Assume the nuclear target mass to be an average between proton and neutron
  constexpr double nucleon_mass = 0.93891875433; // [GeV]
  constexpr double nucleon_mass2 = nucleon_mass * nucleon_mass;

  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if (i == 0) {
      // if proton, skip
      xF[i] = -9999.;
      continue;
    }
    else if (ancestors[i - 1].nucleus == 0 || ancestors[i].nucleus == 1000000000) {
      // if decay process, skip
      // The two conditions should cover differences between Geant4-9 and Geant4-10
      xF[i] = -9999.;
      continue;
    }
    // else if (ancestors[i].pdg == ancestors[i-1].pdg) {
    //   // elastic scatter?
    //   xF[i] = -9999.;
    //   continue;
    // }
    // auto const &inc = ancestors[i - 1];
    auto const& prod = ancestors[i];

    const double mass_inc = ancestor_masses[i - 1];
    const double mass_inc2 = mass_inc * mass_inc;
    const double mass_prod = ancestor_masses[i];

    auto const p_prod = std::sqrt(prod.startpx * prod.startpx + prod.startpy * prod.startpy +
                                  prod.startpz * prod.startpz); // [GeV/c]
    auto const p_inc = std::sqrt(prod.pprodpx * prod.pprodpx + prod.pprodpy * prod.pprodpy +
                                 prod.pprodpz * prod.pprodpz); // [GeV/c]

    double costh =
      (prod.startpx * prod.pprodpx + prod.startpy * prod.pprodpy + prod.startpz * prod.pprodpz) /
      (p_prod * p_inc);

    if (costh > 1.)
      costh = 1.;
    else if (costh < -1.)
      costh = -1.;

    // Calculate the produced particle's longitudinal momentum in the lab frame
    const double pz = p_prod * costh;

    // Calculate the incident particle's energy in the lab frame
    const double E_lab_inc = std::sqrt(p_inc * p_inc + mass_inc2);

    // The center of mass energy of the incident particle and the target nucleon
    const double Ecm = std::sqrt(mass_inc2 + nucleon_mass2 + 2.0 * nucleon_mass * E_lab_inc);
    const double betacm = std::sqrt(E_lab_inc * E_lab_inc - mass_inc2) / (E_lab_inc + nucleon_mass);
    const double gammacm = 1.0 / std::sqrt(1.0 - betacm * betacm); // Lorentz factor

    const double E_lab_prod = std::sqrt(p_prod * p_prod + mass_prod * mass_prod);

    // Boost the produced particle's longitudinal momentum to the center of mass frame
    const double pL = gammacm * (pz - betacm * E_lab_prod);

    // Calculate the Feynman-X (ratio of the produced particle's longitudinal
    // momentum to half of the center of mass energy; the theoretical maximum.)
    xF[i] = 2.0 * pL / Ecm;
  }

  return xF;
}

double calc_weight(const bsim::Decay& decay,
                   const double det_angle,
                   const double energy_ratio,
                   const double parent_energy,
                   const double nu_energy,
                   const double gamma,
                   const ROOT::RVec<double>& rr)
{
  double wght = (det_angle * decay.nimpwt * (energy_ratio * energy_ratio));

  if (std::abs(decay.ptype) != 13) { return wght; }

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
  partial = gamma * (beta[0] * decay.muparpx + beta[1] * decay.muparpy + beta[2] * decay.muparpz);
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

    if (costh > 1.)
      costh = 1.;
    else if (costh < -1.)
      costh = -1.;

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
