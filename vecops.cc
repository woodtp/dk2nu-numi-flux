#include "vecops.h"

double vecops::mag(const std::array<double, 3>& u)
{
  return std::sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]));
}

double vecops::dot(const std::array<double, 3>& u, const std::array<double, 3>& v)
{
  return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

std::array<double, 3> vecops::cross(const std::array<double, 3>& u, const std::array<double, 3>& v)
{
  const std::array<double, 3> res = {
    u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  return res;
}

std::array<double, 3> vecops::project(const std::array<double, 3>& u, const std::array<double, 3>& v)
{
  const double vMag = vecops::mag(v);
  const double costh = vecops::dot(u, v) / (vMag * vMag);
  const std::array<double, 3> res = {u[0] - costh * v[0], u[1] - costh * v[1], u[2] - costh * v[2]};
  return res;
}

double vecops::calculateTheta(const std::array<double, 3>& u, const std::array<double, 3>& v)
{
  const double uMag = vecops::mag(u);
  const double vMag = vecops::mag(v);
  const double costh = vecops::dot(u, v) / (uMag * vMag);
  if (costh > 1.) { return std::acos(1.); }
  else if (costh < -1.) {
    return std::acos(-1.);
  }
  return std::acos(costh); // [rad]
}
