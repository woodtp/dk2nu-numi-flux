#include "dk2nu/tree/dk2nu.h"
#include "ROOT/RVec.hxx"
#include <iostream>


std::vector<std::string> get_volumes(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> get_incident_momenta(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> get_produced_momenta(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> calc_theta(const ROOT::RVec<bsim::Ancestor>&);
std::vector<double> calc_pT(const ROOT::RVec<bsim::Ancestor>&);
std::vector<double> calc_xF(const ROOT::RVec<bsim::Ancestor>&, const ROOT::RVec<double>&);
double calc_weight(const bsim::Decay& decay,
                   double det_angle,
                   double energy_ratio,
                   double parent_energy,
                   double nu_energy,
                   double gamma,
                   const ROOT::RVec<double>& rr);
