import os
import logging

import ROOT
import numpy as np


logging.info("jit'ing functions...")


def set_ROOT_opts(debug: bool = False) -> None:
    logging.info("Setting ROOT options and loading libraries")
    ROOT.gROOT.SetBatch(True)
    if not debug:
        ROOT.EnableImplicitMT()
    dk2nu_lib = os.environ["DK2NU_LIB"]
    ROOT.gSystem.Load(f"{dk2nu_lib}/libdk2nuTree.so")
    ROOT.gSystem.Load("./libWeight.so")
    ROOT.gInterpreter.Declare('#include "Weight.h"')


@ROOT.Numba.Declare(["int"], "double")
def pdg_to_mass(pdg: int) -> float:
    """
    Get the mass of a particle given its PDG code.

    Parameters
    ----------
    pdg : int
        PDG code of the particle.

    Returns
    -------
    mass : float
        Mass of the particle in GeV.
    """

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
    if abs(pdg) == 2212:
        return 0.938272
    if abs(pdg) == 2112:
        return 0.939565
    if abs(pdg) == 3122:
        return 1.115683
    if abs(pdg) == 3222:
        return 1.18937
    if abs(pdg) == 3312:
        return 1.32171
    return -1.0

@ROOT.Numba.Declare(["RVec<double>", "RVec<double>", "RVec<double>"], "RVec<double>")
def calc_magnitudes(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> np.ndarray:
    return np.sqrt(vx**2 + vy**2 + vz**2)


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


@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")
def parent_to_code(parents: np.ndarray) -> np.ndarray:
    def codes(pdg):
        if pdg == 2212:
            return 0
        if pdg == 2112:
            return 1
        if pdg == 211:
            return 2
        if pdg == -211:
            return 3
        if pdg == 321:
            return 4
        if pdg == -321:
            return 5
        if pdg == 130:
            return 6
        if pdg == 13:
            return 7
        if pdg == -13:
            return 8
        else:
            return 9

    return np.array([codes(p) for p in parents[:-1]], dtype=np.int32)


@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")
def target_to_code(targets: np.ndarray) -> np.ndarray:
    def _codes(pdg):
        if pdg == 1000060120:  # carbon
            return 0
        if pdg == 1000130270:  # aluminum
            return 1
        if pdg == 1000260560:  # iron
            return 2
        if pdg == 0 or pdg == 1000000000:  # start process or decay
            return 4
        else:
            return 3

    # account for the apparent Geant4 convention change
    if targets[0] == 0:
        return np.array([_codes(p) for p in targets[1:]], dtype=np.int32)
    return np.array([_codes(p) for p in targets[:-1]], dtype=np.int32)


@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")
def ancestor_parent_pdg(ancestor_pdg: np.ndarray) -> np.ndarray:
    return np.append(np.array([0], dtype=np.int32), ancestor_pdg[:-2])


@ROOT.Numba.Declare(["RVec<int>"], "RVec<double>")
def ancestor_pdg2mass(pdg_ids: np.ndarray):
    """I hate this but it works."""

    def _pdg_to_mass(pdg):
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
        if abs(pdg) == 2212:
            return 0.938272
        if abs(pdg) == 2112:
            return 0.939565
        if abs(pdg) == 3122:
            return 1.115683
        if abs(pdg) == 3222:
            return 1.18937
        if abs(pdg) == 3212:
            return 1.92642
        if abs(pdg) == 3112:
            return 1.19745
        if abs(pdg) == 3312:
            return 1.32171
        if abs(pdg) == 331:
            return 0.95778
        if abs(pdg) == 221:
            return 0.547862
        return -1.0

    return np.array([_pdg_to_mass(p) for p in pdg_ids], dtype=np.float64)
