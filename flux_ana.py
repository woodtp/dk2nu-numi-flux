#!/usr/bin/env python3

import sys
import logging
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

import ROOT  # type: ignore

from root_declarations import set_ROOT_opts
from spectra_definitions import apply_defs


FILE_SETS = {
    "nominal_fhc":  "/pnfs/icarus/persistent/users/awood/g4numi_checks/FHC/false/g4numi*.root",
    "nominal_rhc":  "/pnfs/icarus/persistent/users/awood/g4numi_checks/RHC/false/g4numi*.root",
    "g3Chase_fhc":  "/pnfs/icarus/persistent/users/awood/g4numi_checks/FHC/true/g4numi*.root",
    "g3Chase_rhc":  "/pnfs/icarus/persistent/users/awood/g4numi_checks/RHC/true/g4numi*.root",
    "g4Update_fhc": "/exp/icarus/data/users/awood/uboone_beamsim_g4.10.4/me000z200i/run0/files/g4numi*.root",
    "g4Update_rhc": "/exp/icarus/data/users/awood/uboone_beamsim_g4.10.4/me000z-200i/run0/files/g4numi*.root",
}

POT_PER_FILE = 500_000

ICARUS = [450.37, 7991.98, 79512.66]
TWOXTWO = [0, 0, 103648.837]


def get_pot(files: str) -> int:
    path, pattern = files.rsplit("/", 1)
    nfiles = len(list(Path(path).glob(pattern)))
    return POT_PER_FILE * nfiles


def run_analysis(in_fname: str, out_fname: str, tree_name: str) -> None:
    df = ROOT.RDataFrame("dk2nuTree", in_fname)

    df = apply_defs(df, ICARUS)

    branches = [
        "nu_pdg",
        "parent_pdg",
        "weight",
        "nu_energy",
        "vx",
        "vy",
        "vz",
        "ppvx",
        "ppvy",
        "ppvz",
        "pdpx",
        "pdpy",
        "pdpz",
        "theta_p",
        "par_codes",
        "target_codes",
    ]

    tree_log_str = f"Preparing Tree '{tree_name}' with branches:\n\n"
    for branch in branches:
        tree_log_str += "* " + branch + "\n"

    logging.info(tree_log_str)

    logging.info(
        f"Snapshotting to {out_fname}. Event loop will be executed now, this might take a while..."
    )

    opts = ROOT.RDF.RSnapshotOptions()
    opts.fMode = "UPDATE"

    df.Snapshot(tree_name, out_fname, branches, opts)  # type: ignore


def load_file(fname: str) -> ROOT.RDataFrame:
    return ROOT.RDataFrame("fluxTree", fname)


def main() -> None:
    out_fname = Path("test.root")

    overwrite = True
    file_exists = out_fname.exists()

    if file_exists and not overwrite:
        logging.info(f"{out_fname} already exists and overwriting disabled. Exiting...")
        sys.exit()
    elif file_exists and overwrite:
        logging.info(f"Overwriting {out_fname}")
        out_fname.unlink()

    set_ROOT_opts()
    for horn, files in FILE_SETS.items():
        tree_name = f"fluxTree_{horn}"
        run_analysis(files, str(out_fname), tree_name)

    for horn, files in FILE_SETS.items():
        pot = get_pot(files)
        hpot = ROOT.TH1D(f"hpot_{horn}", "POT", 1, 0, 1)
        hpot.SetBinContent(1, pot)
        logging.info(f"POT = {pot:_}. Writing to {out_fname}:{hpot.GetName()}")
        with ROOT.TFile.Open(str(out_fname), "UPDATE") as _:
            hpot.Write()


if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"\nFinished in {end - start:0.2f} s")
