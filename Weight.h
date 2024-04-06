#include "dk2nu/tree/dk2nu.h"
#include "ROOT/RVec.hxx"
#include <iostream>


double calc_weight(const bsim::Decay& decay,
              const double det_angle,
              const double energy_ratio,
              const double parent_energy,
              const double nu_energy,
              const double gamma,
              const ROOT::RVec<double>& rr
              );
