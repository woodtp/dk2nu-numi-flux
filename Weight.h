#include "dk2nu/tree/dk2nu.h"
#include "ROOT/RVec.hxx"
#include <iostream>

std::vector<std::string> get_volumes(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<bool> is_carbon_vol(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<std::string> interaction_volumes(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> get_incident_momenta(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> get_produced_momenta(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> calc_theta(const ROOT::RVec<bsim::Ancestor>&);
std::vector<double> calc_pT(const ROOT::RVec<bsim::Ancestor>&);
std::vector<double> calc_xF(const ROOT::RVec<bsim::Ancestor>&, const ROOT::RVec<double>&);
double neutrino_energy(const bsim::Decay& decay, const ROOT::RVec<double>& pos);
double calc_weight(const bsim::Decay& decay, const ROOT::RVec<double>& pos);
ROOT::RVec<bool> is_carbon_volume(const ROOT::RVec<bsim::Ancestor>& ancestors);
ROOT::RVec<double> to_rvec(const std::vector<double>& vec);
ROOT::RVec<double> calc_xF_rvec(const ROOT::RVec<bsim::Ancestor>& ancestors,
                                const ROOT::RVec<double>& ancestor_masses);
ROOT::RVec<double> calc_pT_rvec(const ROOT::RVec<bsim::Ancestor>& ancestors);
ROOT::RVec<double> calc_theta_rvec(const ROOT::RVec<bsim::Ancestor>& ancestors);
ROOT::RVec<double> get_incident_momenta_rvec(const ROOT::RVec<bsim::Ancestor>& ancestors);
ROOT::RVec<double> get_produced_momenta_rvec(const ROOT::RVec<bsim::Ancestor>& ancestors);
