#include "dk2nu/tree/dk2nu.h"
#include "ROOT/RVec.hxx"
#include <iostream>


std::vector<std::string> get_volumes(const ROOT::RVec<bsim::Ancestor>& ancestors);
std::vector<double> calc_pT(const ROOT::RVec<bsim::Ancestor>&);
std::vector<double> calc_xF(const ROOT::RVec<bsim::Ancestor>&, const ROOT::RVec<double>&);
double calc_weight(const bsim::Decay& decay,
              const double det_angle,
              const double energy_ratio,
              const double parent_energy,
              const double nu_energy,
              const double gamma,
              const ROOT::RVec<double>& rr
              );
