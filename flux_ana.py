#!/usr/bin/env python3

import argparse
import datetime
import logging
import sys
import time
from pathlib import Path
from typing import Any

import ROOT  # type: ignore
import toml

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

# from spectra_definitions import apply_defs


def get_pot(files: str, pot_per_file: int) -> int:
    path, pattern = files.rsplit("/", 1)
    nfiles = len(list(Path(path).glob(pattern)))
    return pot_per_file * nfiles


def run_analysis(
    in_fname: str,
    out_fname: str,
    tree_name: str,
    cfg: dict[str, Any],
    debug: bool = False,
) -> None:
    df = ROOT.RDataFrame("dk2nuTree", in_fname)  # type: ignore
    if debug:
        df = df.Range(0, 10000)

    root_version = ROOT.__version__  # type: ignore
    minor = root_version.split(".")[1]  # major, minor, patch

    if "/" in minor:
        minor = minor.split("/")[0]

    if int(minor) >= 30:
        ROOT.RDF.Experimental.AddProgressBar(df)  # type: ignore

    logging.debug(f"Loaded {in_fname}. Applying definitions...")

    det_loc: list[float] = [float(x) for x in cfg["location"]]
    if len(det_loc) != 3:
        logging.error(
            f"Invalid location format: {det_loc}. Please provide '[x, y, z]'. Exiting..."
        )
        sys.exit(1)
    logging.info(f"Going to calculate weights for location: {det_loc}")

    logging.debug(f"Location: {det_loc}")

    df = df.Define("det_loc", f"ROOT::RVec<double>{{ {det_loc[0]}, {det_loc[1]}, {det_loc[2]} }}")

    for key, val in cfg["aliases"].items():
        logging.debug(f"Applying Alias: {val} -> {key}")
        df = df.Alias(key, val)

    for key, val in cfg["definitions"].items():
        logging.debug(f"Applying definition: {key}")
        df = df.Define(key, val)

    for val in cfg["filters"].values():
        logging.debug(f"Applying filter: {val}")
        df = df.Filter(val)

    # opts.fCompressionAlgorithm = ROOT.kLZMA
    # opts.fCompressionLevel = 9

    # TODO: Allow for histograms to be defined in the config file.
    # df = df.Filter("nu_pdg == 12")
    #
    # hists = [
    #     df.Histo1D(("hnom_numu", "#nu_{e}", 60, 0, 6), "nu_energy", "weight"),
    #     # df.Filter("parent_pdg == 211").Histo1D(("hnom_numu_pip", "#pi^{+} #to #nu_{e}", 60, 0, 6), "nu_energy", "weight"),
    #     df.Filter("parent_pdg == 321").Histo1D(("hnom_numu_Kp", "K^{+} #to #nu_{e}", 60, 0, 6), "nu_energy", "weight"),
    #     df.Filter("parent_pdg == -13").Histo1D(("hnom_numu_mup", "#mu^{+} #to #nu_{e}", 60, 0, 6), "nu_energy", "weight"),
    # ]
    #
    # with ROOT.TFile.Open(out_fname, "UPDATE") as _:  # type: ignore
    #     for h in hists:
    #         logging.info(f"Writing histogram {h.GetName()} to {out_fname}...")
    #         h.Write()

    tree_log_str = f"Preparing Tree '{tree_name}' with {len(cfg['save_branches'])} branches:\n\n"
    for branch in cfg["save_branches"]:
        tree_log_str += "* " + branch + "\n"

    opts = ROOT.RDF.RSnapshotOptions()  # type: ignore
    opts.fMode = "UPDATE"

    logging.info(tree_log_str)

    logging.info(
        f"Snapshotting to {out_fname}. Event loop will be executed now, this might take a while..."
    )

    if len(cfg["save_branches"]) == 1 and cfg["save_branches"][0] == "*":
        df.Snapshot(tree_name, out_fname, list(cfg["aliases"].keys()) + list(cfg["definitions"].keys()), opts)
    elif len(cfg["save_branches"]) > 0:
        df.Snapshot(tree_name, out_fname, cfg["save_branches"], opts)  # type: ignore
    else:
        logging.info("Not saving any branches. Exiting...")


def load_file(fname: str) -> ROOT.RDataFrame:  # type: ignore
    return ROOT.RDataFrame("fluxTree", fname)  # type: ignore


def main() -> None:
    desc = "Reads Dk2Nu format, calculates weights for the input position, and writes a tree to a new ROOT file."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-c",
        "--config",
        help="Path to configuration file. If unspecified, defaults to the `config.toml` in the current directory.",
        default="config.toml",
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        help="Overwrite output file if it exists",
        action="store_true",
    )
    parser.add_argument("--mt", help="Use multithreading", action="store_true")
    parser.add_argument("-d", "--debug", action="store_true", help="run in debug mode")

    args = parser.parse_args()

    if args.debug:
        logger = logging.getLogger()
        logger.level = logging.DEBUG

    cfg_file = Path(args.config)
    if not cfg_file.exists():
        logging.error("Could not find config.toml. Exiting...")
        sys.exit(1)

    cfg = toml.load(args.config)

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

    if args.debug and args.mt:
        logging.warning(
            "Debug mode and multithreading are not compatible. Disabling multithreading."
        )
        args.mt = False

    set_ROOT_opts(args.mt)

    for name, files in cfg["file_sets"].items():
        parent_dir, pattern = files.rsplit("/", 1)
        n_files = len(list(Path(parent_dir).glob(pattern)))
        logging.info(f"Found {n_files} files in {parent_dir}")
        if n_files == 0:
            logging.error(f"Skipping {name} as no files were found...")
            continue

        tree_name = f"fluxTree_{name}"
        run_analysis(
            files, str(out_fname), tree_name, cfg, debug=args.debug
        )

        pot = get_pot(files, cfg["pot_per_file"])
        hpot = ROOT.TH1D(f"hpot_{name}", "POT", 1, 0, 1)  # type: ignore
        hpot.SetBinContent(1, pot)
        logging.info(f"POT = {pot:_}. Writing to {out_fname}:{hpot.GetName()}")
        with ROOT.TFile.Open(str(out_fname), "UPDATE") as _:  # type: ignore
            hpot.Write()


if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"\nFinished in {datetime.timedelta(seconds=end - start)}")
