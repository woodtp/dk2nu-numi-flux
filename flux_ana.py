#!/u/bin/env python3

import time
from pathlib import Path

import ROOT
# import numpy as np
# import mplhep as hep
# import matplotlib.pyplot as plt

from root_declarations import set_ROOT_opts
from spectra_definitions import apply_defs, build_histograms


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


# def run_analysis(df: ROOT.RDataFrame, det_loc: list[float]) -> ROOT.RDataFrame:
#     rvec_init_str = f"ROOT::RVecD{{ {det_loc[0]}, {det_loc[1]}, {det_loc[2]} }}"
#
#     df = (
#         df.Define("det_loc", rvec_init_str)
#         .Define(
#             "decay_vertex",
#             "ROOT::RVecD{dk2nu.decay.vx, dk2nu.decay.vy, dk2nu.decay.vz}",
#         )
#         .Define("rr", "det_loc - decay_vertex")
#         .Define(
#             "parent_momentum",
#             "ROOT::RVecD{dk2nu.decay.pdpx, dk2nu.decay.pdpy, dk2nu.decay.pdpz}",
#         )
#         .Define("sangdet", "Numba::calc_solid_angle(rr)")
#         .Define("costh", "Numba::calc_costheta_par(parent_momentum, rr)")
#         .Define("parent_mass", "Numba::pdg_to_mass(dk2nu.decay.ptype)")
#         .Define(
#             "parent_energy",
#             "Numba::calc_energy(parent_momentum, parent_mass)",
#         )
#         .Define("pgamma", "Numba::calc_gamma(parent_energy, parent_mass)")
#         .Define("emrat", "Numba::calc_energy_in_beam(pgamma, costh)")
#         .Define("nu_energy", "dk2nu.decay.necm * emrat")
#         .Define(
#             "weight",
#             "calc_weight(dk2nu.decay, sangdet, emrat, parent_energy, nu_energy, pgamma, rr)",
#         )
#         .Define("theta_p", "Numba::theta_p(parent_momentum, det_loc)")
#         .Define("par_codes", "Numba::parent_to_code(dk2nu.ancestor.pdg)")
#         .Define("target_codes", "Numba::target_to_code(dk2nu.ancestor.nucleus)")
#     )
#     return df
#

def main() -> None:
    set_ROOT_opts(ICARUS)
    glob = "g4numi*.root"
    pot = get_pot(UPDATED_G4_FILES["fhc"], glob)

    df = ROOT.RDataFrame(
        "dk2nuTree",
        f"{UPDATED_G4_FILES['fhc']}/{glob}",
    )

    df = apply_defs(df, ICARUS)
    spectra = build_histograms("fhc", "numu")

    print(f"POT: {pot:_}")

    numu = df.Filter("(dk2nu.decay.ntype == 14) && (nu_energy > 0.4)")
    numu_count = numu.Sum("weight").GetValue()
    print(f"numu count: {numu_count}")

    H = []
    for name, hist in spectra.items():
        H.append((name, hist(df)))

    hpot = ROOT.TH1D("hpot", "POT", 1, 0, 1)
    hpot.SetBinContent(1, pot)
    with ROOT.TFile.Open("test.root", "RECREATE") as f:
        for name, h in H:
            h.Draw()  # make sure it's ready
            f.WriteObject(h, name)
        f.WriteObject(hpot, "hpot")


    # h = numu.Histo2D(
    #     ("ints", "ints", 5, 0, 5, 9, 0, 9), "par_codes", "target_codes", "weight"
    # )
    #
    # # print(h_tot)
    #
    # h.Scale(1.0 / numu_count)
    #
    # hep.hist2dplot(h, flow=None)
    # plt.savefig("test.pdf")


if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"\nFinished in {end - start:0.2f} s")
