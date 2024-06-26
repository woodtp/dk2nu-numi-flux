#!/usr/bin/env python3

import argparse
import datetime
import logging
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

import ROOT  # type: ignore
import toml

from spectra_definitions import apply_defs


def get_pot(files: str, pot_per_file: int) -> int:
    path, pattern = files.rsplit("/", 1)
    nfiles = len(list(Path(path).glob(pattern)))
    return pot_per_file * nfiles


def run_analysis(
    in_fname: str,
    out_fname: str,
    tree_name: str,
    location: list[float],
    debug: bool = False,
) -> None:
    df = ROOT.RDataFrame("dk2nuTree", in_fname)
    if debug:
        df = df.Range(0, 1000)

    logging.debug(f"Loaded {in_fname}. Applying definitions...")
    logging.debug(f"Location: {location}")
    df = apply_defs(df, location)

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
        "parent_momentum",
        "theta_p",
        "par_codes",
        "target_codes",
        "ancestor_parent_pdg",
        "ancestor_pT",
        "ancestor_xF",
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
    desc = "Reads Dk2Nu format, calculates weights for the input position, and writes a tree to a new ROOT file."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-c", "--config", help="path to configuration file")
    parser.add_argument(
        "-f",
        "--overwrite",
        help="overwrite output file if it exists",
        action="store_true",
    )
    parser.add_argument("-d", "--debug", action="store_true", help="run in debug mode")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.debug:
        logger = logging.getLogger()
        logger.level = logging.DEBUG

    cfg = toml.load(args.config)

    loc = cfg["location"]
    if len(loc) != 3:
        logging.error(f"Invalid location format: {loc}. Please provide '[x, y, z]'. Exiting...")
        sys.exit(1)

    from root_declarations import set_ROOT_opts

    out_fname = Path(cfg["output_file"])

    file_exists = out_fname.exists()

    if file_exists and not args.overwrite:
        logging.info(
            f"{out_fname} already exists and overwriting disabled. Specify the '-f' flag if you wish to overwite. Exiting..."
        )
        sys.exit()
    elif file_exists and args.overwrite:
        logging.info(f"Overwriting {out_fname}")
        out_fname.unlink()

    set_ROOT_opts(args.debug)

    for name, files in cfg["file_sets"].items():
        tree_name = f"fluxTree_{name}"
        run_analysis(
            files, str(out_fname), tree_name, cfg["location"], debug=args.debug
        )

    for name, files in cfg["file_sets"].items():
        pot = get_pot(files, cfg["pot_per_file"])
        hpot = ROOT.TH1D(f"hpot_{name}", "POT", 1, 0, 1)
        hpot.SetBinContent(1, pot)
        logging.info(f"POT = {pot:_}. Writing to {out_fname}:{hpot.GetName()}")
        with ROOT.TFile.Open(str(out_fname), "UPDATE") as _:
            hpot.Write()


if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"\nFinished in {datetime.timedelta(seconds=end - start)}")
