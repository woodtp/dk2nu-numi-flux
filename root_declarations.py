import logging
import os
import sys
from pathlib import Path

import numpy as np
import ROOT  # type: ignore

logging.info("jit'ing functions...")

def set_ROOT_opts(mt: bool = False) -> None:
    logging.info("Setting ROOT options and loading libraries")
    if mt:
        logging.warning("Enabling ROOT's implicit multithreading. Sometimes this causes issues...")
        ROOT.EnableImplicitMT(10)  # type: ignore

    # libdk2nuTree = Path(__file__).parent / "dk2nu/build/lib/libdk2nuTree.so"
    # if not libdk2nuTree.exists():
    #     logging.error(f"Could not find libdk2nuTree.so at {libdk2nuTree}")
    #     sys.exit(1)

    # ROOT.gSystem.Load(str(libdk2nuTree))  # type: ignore
    dk2nu_libdir = os.environ.get("DK2NU_LIB")
    # logging.debug(dk2nu_dir)
    ROOT.gSystem.Load(f"{dk2nu_libdir}/libdk2nuTree.so")  # type: ignore
    # ROOT.gSystem.AddIncludePath(f" -I. -I{os.environ.get('DK2NU_INC')}")
    # ROOT.gSystem.AddLinkedLibs(f" -L. -L{dk2nu_libdir} -ldk2nuTree -lWeight")
    ROOT.gSystem.Load("libWeight.so")  # type: ignore
    ROOT.gInterpreter.Declare('#include "Weight.h"')  # type: ignore
    ROOT.TH1.AddDirectory(False)  # type: ignore


# @ROOT.Numba.Declare(["int"], "double")  # type: ignore
# def pdg_to_mass(pdg: int) -> float:
#     """
#     Get the mass of a particle given its PDG code.
#
#     Parameters
#     ----------
#     pdg : int
#         PDG code of the particle.
#
#     Returns
#     -------
#     mass : float
#         Mass of the particle in GeV.
#     """
#
#     if abs(pdg) == 13:
#         return 0.1056583755  # muon mass
#     if abs(pdg) == 15:
#         return 1.77686  # tau mass
#     if abs(pdg) == 211:
#         return 0.13957061  # pion mass
#     if abs(pdg) == 321:
#         return 0.493677  # charged kaon mass
#     if pdg == 311 or pdg == 130 or pdg == 310:
#         return 0.497611  # K0 mass
#     if abs(pdg) == 2212:
#         return 0.938272
#     if abs(pdg) == 2112:
#         return 0.939565
#     if abs(pdg) == 3122:
#         return 1.115683
#     if abs(pdg) == 3222:
#         return 1.18937
#     if abs(pdg) == 3312:
#         return 1.32171
#     return -1.0


@ROOT.Numba.Declare(["RVec<double>"], "double")  # type: ignore
def calc_magnitude(v: np.ndarray) -> float:
    """
    Calculate the magnitude of a 3-vector.

    Parameters
    ----------
    v : array_like
        3-vector.

    Returns
    -------
    magnitude : float
        Magnitude of the vector.
    """
    return np.sqrt(np.sum(v**2))

@ROOT.Numba.Declare(["RVec<double>", "RVec<double>", "RVec<double>"], "RVec<double>")  # type: ignore
def calc_magnitudes(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> np.ndarray:
    """
    Calculate the magnitudes given arrays of 3-vector components.

    Parameters
    ----------
    vx : array_like
        array of x-components
    vy : array_like
        array of y-components
    vz : array_like
        array of z-components

    Returns
    -------
    magnitude : array_like
        Magnitudes of the vectors.
    """
    return np.sqrt(vx**2 + vy**2 + vz**2)


# @ROOT.Numba.Declare(["RVec<double>", "double"], "double")  # type: ignore
# def calc_energy(p: np.ndarray, mass: float) -> float:
#     """Calculate the energy of a particle given from its four-momentum."""
#     return np.sqrt(np.sum(p**2) + mass**2)


# @ROOT.Numba.Declare(["double", "double"], "double")  # type: ignore
# def calc_gamma(energy: float, mass: float) -> float:
#     """Calculate the gamma factor of a particle given its energy and mass."""
#     return energy / mass


# @ROOT.Numba.Declare(["double"], "double")  # type: ignore
# def calc_beta(gamma: float) -> float:
#     return np.sqrt(1.0 - 1.0 / gamma**2)


# @ROOT.Numba.Declare(["double", "double"], "double")  # type: ignore
# def calc_energy_in_beam(gamma: float, costh: float) -> float:
#     """ Calculate the fraction of the energy of a particle in the beam direction.
#     """
#     beta = np.sqrt(1 - 1 / gamma**2)
#     return 1 / (gamma * (1 - beta * costh))


# @ROOT.Numba.Declare(["RVec<double>"], "double")  # type: ignore
# def calc_solid_angle(rr: np.ndarray) -> float:
#     """
#     Calculate solid angle/4pi of a point (small angle approximation).
#
#     Parameters
#     ----------
#     rr : array_like
#         3-vector of the position vector of the point.
#         Typically the vector from the decay vertex to the detector location.
#
#     Returns
#     -------
#     solid_angle : float
#         solid angle of the point.
#     """
#     rdet2 = 100**2  # for converting to cm
#     rrmag2 = rr[0]**2 + rr[1]**2 + rr[2]**2
#     return rdet2 / rrmag2 / (4.0 * np.pi)


# @ROOT.Numba.Declare(["RVec<double>", "RVec<double>"], "double")  # type: ignore
# def calc_costheta_par(parent_momentum: np.ndarray, rr: np.ndarray) -> float:
#     """
#     Calculate the angle between a momentum vector and position. If either the parent momentum or the position vector is
#     zero, the function returns -9999.0. Otherwise, the result is clipped to [-1, 1].
#
#     Parameters
#     ----------
#     parent_momentum : array_like
#         3-vector of the parent particle's momentum.
#     rr : array_like
#         3-vector of the position vector of the point.
#
#     Returns
#     -------
#     costh : float
#         boost factor
#     """
#
#     if not parent_momentum.any() or not rr.any():
#         return -9999.0
#
#     mom_mag = np.sqrt(np.sum(parent_momentum**2))
#     rr_mag = np.sqrt(np.sum(rr**2))
#     dot = parent_momentum[0]*rr[0] + parent_momentum[1]*rr[1] + parent_momentum[2]*rr[2]
#
#     return max(min(dot / (mom_mag * rr_mag), 1.0), -1.0)



@ROOT.Numba.Declare(["RVec<double>", "RVec<double>"], "double")  # type: ignore
def theta_p(p: np.ndarray, det_loc: np.ndarray) -> float:
    """
    Calculate the angle (in mrad) between a momentum vector and NuMI within the plane
    connecting the detector location and the NuMI axis.
    If the provided location is colinear with the beam axis, the angle is calculated with respect to the x-z plane.
    If the input momentum vector is zero, the function returns -9999.0.

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
    if not p.any():
        return -9999.0

    numi_axis = np.array([0, 0, 1.0])

    n = np.cross(det_loc, numi_axis)

    if not n.any():
        n = np.array([0, 1.0, 0])

    # computing dot product manually to suppress warnings about
    # the input RVecD not being in contiguous memory.
    proj_mag = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / np.sum(n**2)

    projection = p - proj_mag * n

    costh = (projection @ numi_axis) / np.sqrt(np.sum(projection**2))
    theta = np.arccos(costh)

    crossprod = np.cross(projection, numi_axis)
    if (crossprod @ n) < 0:
        theta = -theta
    return 1000*theta

@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")  # type: ignore
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
            return 8
        if pdg == -13:
            return 9
        else:
            return 7

    return np.array([codes(p) for p in parents[:-1]], dtype=np.int32)


@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")  # type: ignore
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


@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")  # type: ignore
def ancestor_parent_pdg(ancestor_pdg: np.ndarray) -> np.ndarray:
    return np.append(np.array([0], dtype=np.int32), ancestor_pdg[:-2])

@ROOT.Numba.Declare(["RVec<double>"], "RVec<double>")  # type: ignore
def ancestor_incident_momenta(momenta: np.ndarray) -> np.ndarray:
    return momenta[1:]

# @ROOT.Numba.Declare(["RVec<int>"], "RVec<double>")  # type: ignore
# def ancestor_pdg2mass(pdg_ids: np.ndarray):
#     """I hate this but it works."""
#
#     def _pdg_to_mass(pdg):
#         if abs(pdg) == 13:
#             return 0.1056583755  # muon mass
#         if abs(pdg) == 15:
#             return 1.77686  # tau mass
#         if abs(pdg) == 211:
#             return 0.13957061  # pion mass
#         if abs(pdg) == 321:
#             return 0.493677  # charged kaon mass
#         if pdg == 311 or pdg == 130 or pdg == 310:
#             return 0.497611  # K0 mass
#         if abs(pdg) == 2212:
#             return 0.938272
#         if abs(pdg) == 2112:
#             return 0.939565
#         if abs(pdg) == 3122:
#             return 1.115683
#         if abs(pdg) == 3222:
#             return 1.18937
#         if abs(pdg) == 3212:
#             return 1.92642
#         if abs(pdg) == 3112:
#             return 1.19745
#         if abs(pdg) == 3312:
#             return 1.32171
#         if abs(pdg) == 331:
#             return 0.95778
#         if abs(pdg) == 221:
#             return 0.547862
#         return -1.0
#
#     return np.array([_pdg_to_mass(p) for p in pdg_ids], dtype=np.float64)

@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")  # type: ignore
def incident_pdg_codes(ancestor_pdg: np.ndarray) -> np.ndarray:
    return ancestor_pdg[:-2]

@ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")  # type: ignore
def produced_pdg_codes(ancestor_pdg: np.ndarray) -> np.ndarray:
    return ancestor_pdg[1:-1]
