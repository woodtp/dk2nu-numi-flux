#define _USE_MATH_DEFINES
#include <cmath>

#include "Weight.h"

static constexpr double DEFAULT_DOUBLE = -9999.0;
static constexpr double MUON_MASS = 0.1056583755; // [GeV/c^2]
static constexpr double TAU_MASS = 1.77686; // [GeV/c^2]
static constexpr double SQM_TO_SQCM = 100.0*100.0;
static constexpr double FOUR_PI = 4.0 * M_PI;

enum PDG {
  ELECTRON = 11,
  ELECTRON_NEUTRINO = 12,
  MUON = 13,
  MUON_NEUTRINO = 14,
  TAU = 15,
  TAU_NEUTRINO = 16,
  PION = 211,
  KAON = 321,
  K0L = 130,
  K0S = 310,
  K0 =  311,
  P = 2212,
  N = 2112,
  LAMBDA = 3122,
  SIGMAP = 3222,
  SIGMA0 = 3212,
  SIGMAM = 3112,
  XIM = 3312,
  ETA = 221,
  ETAPRIME = 331
};

static double clip(const double val)
{
  // clip value to [-1, 1]
  return std::max(-1.0, std::min(1.0, val));
}

static bool skip_particle(const bsim::Ancestor& ancestor)
{
  // if the beam proton or decay process, skip
  // this should cover differences between Geant4-9 (Primary) and Geant4-10 (BeamParticle)
  return ancestor.proc == "Primary" || ancestor.proc == "BeamParticle" || ancestor.proc == "Decay";
}

std::vector<std::string> get_volumes(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<std::string> vol;
  vol.reserve(ancestors.size() - 1);
  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    vol.push_back(ancestors[i].ivol);
  }
  return vol;
}

std::vector<bool> is_carbon_vol(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<bool> vol;
  vol.reserve(ancestors.size() - 2);
  for (std::size_t i = 1; i < ancestors.size() - 1; ++i) {
    vol.emplace_back(ancestors[i].ivol == "TGT1" || ancestors[i].ivol == "Budal_Monitor" || ancestors[i].ivol == "Budal_HFVS" || ancestors[i].ivol == "Budal_VFHS");
  }
  return vol;
}

std::vector<double> get_incident_momenta(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<double> p_inc(ancestors.size() - 1, DEFAULT_DOUBLE);
  const bool is_old_g4 = ancestors[0].proc == "Primary";
  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if(skip_particle(ancestors[i])) continue;
    const bsim::Ancestor& part = i == 0 || !is_old_g4 ? ancestors[i] : ancestors[i - 1];
    p_inc[i] = std::sqrt(part.pprodpx*part.pprodpx +
                         part.pprodpy*part.pprodpy +
                         part.pprodpz*part.pprodpz);
  }
  return p_inc;
}

std::vector<double> get_produced_momenta(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<double> p_prod(ancestors.size() - 1, DEFAULT_DOUBLE);
  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if(skip_particle(ancestors[i])) continue;
    p_prod[i] = std::sqrt(ancestors[i].startpx*ancestors[i].startpx
                          + ancestors[i].startpy*ancestors[i].startpy
                          + ancestors[i].startpz*ancestors[i].startpz);
  }
  return p_prod;
}

std::vector<double> calc_theta(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<double> thetas(ancestors.size() - 1, DEFAULT_DOUBLE);

  const bool is_old_g4 = ancestors[0].proc == "Primary";

  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if (skip_particle(ancestors[i])) continue;

    auto const& prod = ancestors[i];

    const double pincpx = is_old_g4 ? ancestors[i - 1].pprodpx : prod.pprodpx;
    const double pincpy = is_old_g4 ? ancestors[i - 1].pprodpy : prod.pprodpy;
    const double pincpz = is_old_g4 ? ancestors[i - 1].pprodpz : prod.pprodpz;
    const double p_inc = std::sqrt(pincpx*pincpx + pincpy*pincpy + pincpz*pincpz); // [GeV/c]

    const double p_prod = std::sqrt(prod.startpx*prod.startpx +
                                    prod.startpy*prod.startpy +
                                    prod.startpz*prod.startpz); // [GeV/c]


    if (p_prod == 0. || p_inc == 0.) continue;

    const double costh = clip((prod.startpx*pincpx + prod.startpy*pincpy + prod.startpz*pincpz) / (p_prod*p_inc));

    thetas[i] = std::acos(costh);
  }

  return thetas;
}

std::vector<double> calc_pT(const ROOT::RVec<bsim::Ancestor>& ancestors)
{
  std::vector<double> pT(ancestors.size() - 1, DEFAULT_DOUBLE);

  const bool is_old_g4 = ancestors[0].proc == "Primary";

  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if (skip_particle(ancestors[i])) continue;

    auto const& prod = ancestors[i];

    const double pincpx = is_old_g4 ? ancestors[i - 1].pprodpx : prod.pprodpx;
    const double pincpy = is_old_g4 ? ancestors[i - 1].pprodpy : prod.pprodpy;
    const double pincpz = is_old_g4 ? ancestors[i - 1].pprodpz : prod.pprodpz;
    const double p_inc = std::sqrt(pincpx*pincpx + pincpy*pincpy + pincpz*pincpz); // [GeV/c]

    const double p_prod = std::sqrt(prod.startpx*prod.startpx +
                                    prod.startpy*prod.startpy +
                                    prod.startpz*prod.startpz);


    if (p_prod == 0. || p_inc == 0.) continue;

    const double costh = clip((prod.startpx*pincpx + prod.startpy*pincpy + prod.startpz*pincpz) / (p_prod*p_inc));

    const double sinth = std::sqrt(1. - costh*costh);

    pT[i] = p_prod*sinth;
  }

  return pT;
}

std::vector<double> calc_xF(const ROOT::RVec<bsim::Ancestor>& ancestors,
                            const ROOT::RVec<double>& ancestor_masses)
{
  std::vector<double> xF(ancestors.size() - 1, DEFAULT_DOUBLE);

  // Assume the nuclear target mass to be an average between proton and neutron
  constexpr double nucleon_mass = 0.93891875433; // [GeV]
  constexpr double nucleon_mass2 = nucleon_mass*nucleon_mass;

  const bool is_old_g4 = ancestors[0].proc == "Primary";

  for (std::size_t i = 0; i < ancestors.size() - 1; ++i) {
    if (skip_particle(ancestors[i])) continue;

    auto const& prod = ancestors[i];

    const double mass_inc = ancestor_masses[i - 1];
    const double mass_inc2 = mass_inc*mass_inc;
    const double mass_prod = ancestor_masses[i];

    const double pincpx = is_old_g4 ? ancestors[i - 1].pprodpx : prod.pprodpx;
    const double pincpy = is_old_g4 ? ancestors[i - 1].pprodpy : prod.pprodpy;
    const double pincpz = is_old_g4 ? ancestors[i - 1].pprodpz : prod.pprodpz;
    const double p_inc = std::sqrt(pincpx*pincpx + pincpy*pincpy + pincpz*pincpz); // [GeV/c]

    auto const p_prod = std::sqrt(prod.startpx*prod.startpx +
                                  prod.startpy*prod.startpy +
                                  prod.startpz*prod.startpz); // [GeV/c]

    if (p_prod == 0. || p_inc == 0.) continue;

    const double costh = clip((prod.startpx*pincpx +
                               prod.startpy*pincpy +
                               prod.startpz*pincpz) / (p_prod*p_inc));

    // Calculate the produced particle's longitudinal momentum in the lab frame
    const double pz = p_prod*costh;

    // Calculate the incident particle's energy in the lab frame
    const double E_lab_inc = std::sqrt(p_inc*p_inc + mass_inc2);

    // The center of mass energy of the incident particle and the target nucleon
    const double Ecm = std::sqrt(mass_inc2 + nucleon_mass2 + 2.0*nucleon_mass*E_lab_inc);
    const double betacm = std::sqrt(E_lab_inc*E_lab_inc - mass_inc2) / (E_lab_inc + nucleon_mass);
    const double gammacm = 1.0 / std::sqrt(1.0 - betacm*betacm); // Lorentz factor

    const double E_lab_prod = std::sqrt(p_prod*p_prod + mass_prod*mass_prod);

    // Boost the produced particle's longitudinal momentum to the center of mass frame
    const double pL = gammacm*(pz - betacm*E_lab_prod);

    // Calculate the Feynman-X (ratio of the produced particle's longitudinal
    // momentum to half of the center of mass energy; the theoretical maximum.)
    xF[i] = 2.0*pL / Ecm;
  }

  return xF;
}

static double pdgid_to_mass(int pdgid) {
  switch(std::abs(pdgid)) {
    case PDG::MUON:
        return MUON_MASS;
    case PDG::TAU:
        return TAU_MASS;
    case PDG::PION:
        return 0.13957061;
    case PDG::KAON:
        return 0.493677;
    case PDG::K0L:
    case PDG::K0S:
    case PDG::K0:
        return 0.497611;
    case PDG::P:
        return 0.938272;
    case PDG::N:
        return 0.939565;
    case PDG::LAMBDA:
        return 1.115683;
    case PDG::SIGMAP:
        return 1.18937;
    case PDG::SIGMA0:
        return 1.92642;
    case PDG::SIGMAM:
        return 1.19745;
    case PDG::XIM:
        return 1.32171;
    case PDG::ETAPRIME:
        return 0.95778;
    case PDG::ETA:
        return 0.547862;
    default:
      return -1.0;
  }
}

static double parent_energy(const bsim::Decay& decay) {
  const double parent_mass = pdgid_to_mass(decay.ptype);
  const double parent_p2 = decay.pdpx*decay.pdpx + decay.pdpy*decay.pdpy + decay.pdpz*decay.pdpz;
  return std::sqrt(parent_p2 + parent_mass*parent_mass);
}

static double emratio(const bsim::Decay& decay, const ROOT::RVec<double>& pos) {
  if (const auto size = pos.size() != 3) {
    std::cout << "[ERROR:Weight.cc:emratio] Input position had length != 3. pos.size() = " << size << std::endl;
    exit(1);
  }
  const double gamma = parent_energy(decay) / pdgid_to_mass(decay.ptype);
  const double beta = std::sqrt(1 - 1 / (gamma*gamma));

  const ROOT::RVec<double> rr {pos[0] - decay.vx, pos[1] - decay.vy, pos[2] - decay.vz};
  const double rrmag2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];

  const double parent_p2 = decay.pdpx*decay.pdpx + decay.pdpy*decay.pdpy + decay.pdpz*decay.pdpz;
  const double costh = clip((decay.pdpx*rr[0] + decay.pdpy*rr[1] + decay.pdpz*rr[2]) / std::sqrt(parent_p2*rrmag2));
  const double emratio = 1 / (gamma * (1 - beta * costh));

  return emratio;
}

double neutrino_energy(const bsim::Decay& decay, const ROOT::RVec<double>& pos) {
  return decay.necm*emratio(decay, pos);
}

static double dot_product(const double* a, const double* b, const std::size_t dim) {
  double result = 0.;
  for (std::size_t i = 0; i < dim; ++i) {
    result += a[i]*b[i];
  }
  return result;
}

double calc_weight(const bsim::Decay& decay,
                   const ROOT::RVec<double>& pos)
{
  // Yoinked from https://github.com/NuSoftHEP/dk2nu/blob/main/tree/calcLocationWeights.cxx#L48
  if (const auto size = pos.size() != 3) {
    std::cout << "[ERROR:Weight.cc:calc_weight] Input position had length != 3. pos.size() = " << size << std::endl;
    exit(1);
  }

  const ROOT::RVec<double> rr {pos[0] - decay.vx, pos[1] - decay.vy, pos[2] - decay.vz};
  const double rrmag2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];

  const double sangdet = SQM_TO_SQCM / rrmag2 / FOUR_PI;

  const double emrat = emratio(decay, pos);
  const double wght = sangdet*decay.nimpwt*emrat*emrat; // nu flux [m^{-2}]

  // if the nu parent is not a muon, we're done.
  if (std::abs(decay.ptype) != PDG::MUON) return wght;

  // For muon parents, we need to apply correction for anisotropic decay

  double p_energy = parent_energy(decay);
  // Boost new neutrino to mu decay cm
  const double beta[] = {
    decay.pdpx / p_energy,
    decay.pdpy / p_energy,
    decay.pdpz / p_energy
  };

  const double rr_mag = std::sqrt(rrmag2);

  const double nu_energy = decay.necm*emrat;

  const double p_nu[] = {
    rr[0]*nu_energy / rr_mag,
    rr[1]*nu_energy / rr_mag,
    rr[2]*nu_energy / rr_mag
  };

  double gamma = parent_energy(decay) / pdgid_to_mass(decay.ptype);
  double partial = gamma*dot_product(beta, p_nu, 3);
  partial = nu_energy - partial / (gamma + 1.);

  // compute 4-momentum of the neutrino in the muon decay center of mass frame
  double p_dcm_nu[4] = {0.};
  for (int i = 0; i < 3; ++i) {
    p_dcm_nu[i] = p_nu[i] - beta[i]*gamma*partial;
    p_dcm_nu[3] += p_dcm_nu[i]*p_dcm_nu[i];
  }
  p_dcm_nu[3] = std::sqrt(p_dcm_nu[3]);

  // Boost parent of mu to mu production cm
  const double gamma_mu = decay.ppenergy / MUON_MASS;

  const double beta_mu[] = {
    decay.ppdxdz*decay.pppz / decay.ppenergy,
    decay.ppdydz*decay.pppz / decay.ppenergy,
    decay.pppz / decay.ppenergy
  };

  const double muparp[] = {
    decay.muparpx,
    decay.muparpy,
    decay.muparpz
  };

  partial = gamma_mu*dot_product(beta_mu, muparp, 3);
  partial = decay.mupare - partial / (gamma_mu + 1.);

  // compute 4-momentum of the parent of the muon in the muon production center of mass frame
  double p_pcm_mp[4] = {0.};
  for (int i = 0; i < 3; ++i) {
    p_pcm_mp[i] = muparp[i] - beta_mu[i]*gamma_mu*partial;
    p_pcm_mp[3] += p_pcm_mp[i]*p_pcm_mp[i];
  }
  p_pcm_mp[3] = std::sqrt(p_pcm_mp[3]);

  double wt_ratio = 1.;
  //have to check p_pcm_mp
  //it can be 0 if mupar..=0. (I guess muons created in target??)
  if (p_pcm_mp[3] != 0.) {
    //calc new decay angle w.r.t. (anti)spin direction
    const double costh = clip(dot_product(p_dcm_nu, p_pcm_mp, 4) / (p_dcm_nu[3]*p_pcm_mp[3]));

    //calc relative weight due to angle difference
    if (std::abs(decay.ntype) == PDG::ELECTRON_NEUTRINO) { wt_ratio = 1. - costh; }
    else if (std::abs(decay.ntype) == PDG::MUON_NEUTRINO) {
      const double xnu = 2.*decay.necm / MUON_MASS;
      wt_ratio = ((3. - 2.*xnu) - (1. - 2.*xnu)*costh) / (3. - 2.*xnu);
    }
    else if (std::abs(decay.ntype) == PDG::TAU_NEUTRINO) {
      const double xnu = 2.*decay.necm / TAU_MASS;
      wt_ratio = ((3. - 2.*xnu) - (1. - 2.*xnu)*costh) / (3. - 2.*xnu);
      std::cout << "[INFO:Weight.cc:calc_weight] calculating weight for tau neutrino; this may not be correct" << std::endl;
    }
    else {
      std::cout << "[INFO:Weight.cc:calc_weight] eventRates:: Bad neutrino type = " << decay.ntype << std::endl;
    }
  }
  return wght*wt_ratio;
}
