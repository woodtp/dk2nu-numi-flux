import os

import ROOT
import numpy as np


def set_ROOT_opts(det_loc: list[float]) -> None:
    dk2nu_lib = os.environ["DK2NU_LIB"]
    ROOT.gSystem.Load(f"{dk2nu_lib}/libdk2nuTree.so")
    ROOT.gROOT.SetBatch(True)
    ROOT.EnableImplicitMT()
    ROOT.gInterpreter.Declare(
        f"static const std::array<double, 3> DETECTOR_LOCATION = {{{det_loc[0]}, {det_loc[1]}, {det_loc[2]}}};"
    )
    ROOT.gInterpreter.Declare('#include "Weight.h"')


@ROOT.Numba.Declare(["int"], "double")
def pdg_to_mass(pdg: int) -> float:
    if abs(pdg) == 13:
        return 0.13957039  # muon mass
    if abs(pdg) == 15:
        return 1.77686  # tau mass
    if abs(pdg) == 211:
        return 0.13957061  # pion mass
    if abs(pdg) == 321:
        return 0.493677  # charged kaon mass
    if pdg == 311 or pdg == 130 or pdg == 310:
        return 0.497611  # K0 mass
    return -1.0


@ROOT.Numba.Declare(["RVec<double>", "double"], "double")
def calc_energy(p: np.ndarray, mass: float) -> float:
    """Calculate the energy of a particle given from its four-momentum."""
    return np.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2 + mass**2)


@ROOT.Numba.Declare(["double", "double"], "double")
def calc_gamma(energy: float, mass: float) -> float:
    """Calculate the gamma factor of a particle given its energy and mass."""
    return energy / mass


@ROOT.Numba.Declare(["double"], "double")
def calc_beta(gamma: float) -> float:
    return np.sqrt(1.0 - 1.0 / gamma**2)


@ROOT.Numba.Declare(["double", "double"], "double")
def calc_energy_in_beam(gamma: float, costh: float) -> float:
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return 1.0 / (gamma * (1.0 - beta * costh))


@ROOT.Numba.Declare(["RVec<double>"], "double")
def calc_solid_angle(rr: np.ndarray) -> float:
    """
    Calculate solid angle/4pi of a point.

    Parameters
    ----------
    rr : array_like
        3-vector of the position vector of the point.
        Typically the vector from the decay vertex to the detector location.

    Returns
    -------
    solid_angle : float
        solid angle of the point.
    """
    rdet = 100.0  # for converting to cm
    mag2 = rr[0] ** 2 + rr[1] ** 2 + rr[2] ** 2
    return (rdet**2) / mag2 / (4.0 * np.pi)


@ROOT.Numba.Declare(["RVec<double>", "RVec<double>"], "double")
def calc_costheta_par(parent_momentum: np.ndarray, rr: np.ndarray) -> float:
    """
    Calculate a boost correction factor from the parent particle's momentum.

    Parameters
    ----------
    parent_momentum : array_like
        3-vector of the parent particle's momentum.
    rr : array_like
        3-vector of the position vector of the point.

    Returns
    -------
    costh : float
        boost factor
    """

    dot = (
        parent_momentum[0] * rr[0]
        + parent_momentum[1] * rr[1]
        + parent_momentum[2] * rr[2]
    )
    mom_mag = np.sqrt(
        parent_momentum[0] ** 2 + parent_momentum[1] ** 2 + parent_momentum[2] ** 2
    )
    rr_mag = np.sqrt(rr[0] ** 2 + rr[1] ** 2 + rr[2] ** 2)
    costh = dot / (mom_mag * rr_mag)
    if costh > 1.0:
        return 1.0
    if costh < -1.0:
        return -1.0
    return costh


@ROOT.Numba.Declare(["RVec<double>", "RVec<double>"], "double")
def theta_p(p: np.ndarray, det_loc: np.ndarray) -> float:
    """
    Calculate the angle (in degrees) between a momentum vector and NuMI within the plane
    connecting the detector location and the NuMI axis.

    Parameters
    ----------
        p : array_like
            3-vector of the momentum vector.
        det_loc : array_like
            3-vector of the detector location.

    Returns
    -------
         theta : float
         The angle between the momentum vector and the NuMI axis in degrees.
    """

    numi_axis = np.array([0, 0, 1.0])

    n = np.cross(det_loc, numi_axis)

    # computing dot product manually to suppress warnings about
    # the input RVecD not being in contiguous memory.
    proj_mag = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / (n**2).sum()

    projection = p - proj_mag * n

    costh = (projection @ numi_axis) / (
        np.linalg.norm(projection) * np.linalg.norm(numi_axis)
    )
    theta = np.degrees(np.arccos(costh))

    crossprod = np.cross(projection, numi_axis)
    if (crossprod @ n) < 0:
        theta = -theta
    return theta