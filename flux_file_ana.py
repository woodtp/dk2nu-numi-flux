#!/usr/bin/env python3

import time
from pathlib import Path
import sys

import ROOT
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt

from root_declarations import set_ROOT_opts


UPDATED_G4_FILES = {
    "fhc": "/exp/icarus/data/users/awood/uboone_beamsim_g4.10.4/me000z200i/run0/files",
    "rhc": "/exp/icarus/data/users/awood/uboone_beamsim_g4.10.4/me000z-200i/run0/files",
}

POT_PER_FILE = 500_000

ICARUS = [450.37, 7991.98, 79512.66]
TWOXTWO = [0, 0, 103648.837]

def get_pot(file_path: str, pattern: str) -> int:
    nfiles = len(list(Path(file_path).glob(pattern)))
    return POT_PER_FILE * nfiles


def run_analysis(df: ROOT.RDataFrame, det_loc: list[float]) -> ROOT.RDataFrame:
    rvec_init_str = f"ROOT::RVecD{{ {det_loc[0]}, {det_loc[1]}, {det_loc[2]} }}"

    df = (
        df.Define("det_loc", rvec_init_str)
        .Define(
            "decay_vertex",
            "ROOT::RVecD{dk2nu.decay.vx, dk2nu.decay.vy, dk2nu.decay.vz}",
        )
        .Define("rr", "det_loc - decay_vertex")
        .Define(
            "parent_momentum",
            "ROOT::RVecD{dk2nu.decay.pdpx, dk2nu.decay.pdpy, dk2nu.decay.pdpz}",
        )
        .Define("det_angle", "Numba::calc_det_angle(rr)")
        .Define("costh", "Numba::calc_costheta(parent_momentum, rr)")
        .Define("parent_mass", "Numba::pdg_to_mass(dk2nu.decay.ptype)")
        .Define(
            "parent_energy",
            "Numba::calc_energy(parent_momentum, parent_mass)",
        )
        .Define("pgamma", "Numba::calc_gamma(parent_energy, parent_mass)")
        .Define("emrat", "Numba::calc_energy_ratio(pgamma, costh)")
        .Define("nu_energy", "dk2nu.decay.necm * emrat")
        .Define("weight", "calc_weight(dk2nu.decay, det_angle, emrat, parent_energy, nu_energy, pgamma, rr)")
        .Define("theta_p", "Numba::theta_p(parent_momentum, det_loc)")
    )
    return df


def main() -> None:
    set_ROOT_opts(ICARUS)
    glob = "g4numi*.root"
    pot = get_pot(UPDATED_G4_FILES["fhc"], glob)

    df = ROOT.RDataFrame(
        "dk2nuTree",
        f"{UPDATED_G4_FILES['fhc']}/{glob}",
    )

    df = run_analysis(df, ICARUS)

    print(pot)

    h_tot = df.Filter("(dk2nu.decay.ntype == 14) && (dk2nu.decay.ptype == 321)").Histo1D(
        ("nu", "nu", 30, -15, 15), "theta_p", "weight"
    )

    h_tot.Scale(1.0 / pot)

    hep.histplot(h_tot, edges=False)
    plt.savefig("test.pdf")


if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"\nFinished in {end - start:0.2f} s")
