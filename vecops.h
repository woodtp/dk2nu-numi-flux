#ifndef VECOPS_H
#define VECOPS_H

#include <array>
#include <cmath>

namespace vecops {

  double mag(const std::array<double, 3>& u);
  double dot(const std::array<double, 3>& u, const std::array<double, 3>& v);
  std::array<double, 3> cross(const std::array<double, 3>& u, const std::array<double, 3>& v);
  std::array<double, 3> project(const std::array<double, 3>& u, const std::array<double, 3>& v);
  double calculateTheta(const std::array<double, 3>& u, const std::array<double, 3>& v);

} // namespace vecops
#endif
