#!/usr/bin/env python3

import argparse
from collections.abc import Iterable
import datetime
import logging
import sys
import time
from pathlib import Path
from typing import Any, Callable, Optional

import ROOT  # type: ignore
import uproot

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

def get_pot(files: str) -> int:
    df = ROOT.RDataFrame("dkmetaTree", files)
    return df.Sum("pots").GetValue()


def run_analysis(
    in_fname: str,
    out_fname: str,
    sample_name: str,
    det_loc: list[float],
    cfg: dict[str, Any],
    debug: bool = False,
) -> None:
    df = ROOT.RDataFrame("dk2nuTree", in_fname)  # type: ignore
    if debug:
        df = df.Range(0, 250_000)

    root_version = ROOT.__version__  # type: ignore
    minor = root_version.split(".")[1]  # major, minor, patch

    if "/" in minor:
        minor = minor.split("/")[0]

    if int(minor) >= 30:
        ROOT.RDF.Experimental.AddProgressBar(df)  # type: ignore

    logging.debug(f"Loaded {in_fname}. Applying definitions...")

    pot = get_pot(in_fname)
    pot_wgt = 1.0 / pot

    logging.info(f"POT = {pot:_}. Weight available as 'pot_wgt' = 1/POT = {pot_wgt:e}")
    df = df.Define("pot_wgt", str(pot_wgt))

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

    det_loc_str = f"{det_loc[0]}_{det_loc[1]}_{det_loc[2]}"

    def _is_regular_bins(bins):
        return len(bins) == 3 and all(isinstance(x, int) or isinstance(x, float) for x in bins)

    def _is_variable_bins(bins):
        return all(isinstance(x, list) and len(x) == 3 for x in bins)

    hists = []
    df_cache = {}
    for key, h in cfg["histograms"].items():
        if (other_key := h.get("copy_from")) is not None:
            logging.debug(f"Copying histogram params from {other_key} to {key}")
            h = cfg["histograms"][other_key].copy()
            # apply any overrides
            h.update(cfg["histograms"][key])

        weight = h.get("weight", "")

        df_ref = df
        if (filt_def := h.get("filter")) is not None:
            if filt_def in df_cache:
                df_ref = df_cache[filt_def]
            else:
                df_ref = df.Filter(filt_def)
                df_cache[filt_def] = df_ref

        if (dir := h.get("dir")) is not None:
            hist_name = f"{sample_name}/{dir}/{key}"
        else:
            hist_name = f"{sample_name}/{key}"

        # checking whether it's a 1D or 2D histogram
        def _hist_info(h):
            info = "Registering histogram:\n"
            for k, v in h.items():
                info += f"    {k}: {v}\n"
            return info

        logging.info(_hist_info(h))
        if "yvar" not in h:
            if _is_regular_bins(h["bins"]):
                nbins, low, high = h["bins"]

                logging.debug(f"Creating 1D histogram {key} with {nbins} bins from {low} to {high}")

                hists.append(df_ref.Histo1D((hist_name, h["title"], nbins, low, high), h["var"], weight))
            elif _is_variable_bins(h["bins"]):
                bins = np.array([])
                for bs in h["bins"]:
                    low, high, step = bs
                    bins = np.append(bins, np.arange(low, high+step, step))

                logging.debug(f"Creating 1D histogram {key} with variable bins:\n{bins}")

                hists.append(df_ref.Histo1D((hist_name, h["title"], len(bins) - 1, bins), h["var"], weight))
            else:
                raise ValueError(f"Invalid binning format for histogram {key}. Exiting...")
        else:
            if _is_regular_bins(h["xbins"]) and _is_regular_bins(h["ybins"]):
                nbinsx, xlow, xhigh = h["xbins"]
                nbinsy, ylow, yhigh = h["ybins"]

                hists.append(df_ref.Histo2D((hist_name, h["title"], nbinsx, xlow, xhigh, nbinsy, ylow, yhigh), h["xvar"], h["yvar"], weight))
            elif _is_variable_bins(h["xbins"]) and _is_regular_bins(h["ybins"]):
                xbins = np.array([])

                for bs in h["xbins"]:
                    low, high, step = bs
                    xbins = np.append(xbins, np.arange(low, high+step, step))

                nbinsy, ylow, yhigh = h["ybins"]

                logging.debug(f"Creating 2D histogram {key} with variable x bins:\n{xbins} and {nbinsy} regular y bins from {ylow} to {yhigh}")

                hists.append(df_ref.Histo2D((hist_name, h["title"], len(xbins)-1, xbins, nbinsy, ylow, yhigh), h["xvar"], h["yvar"], weight))
            elif _is_regular_bins(h["xbins"]) and _is_variable_bins(h["ybins"]):
                nbinsx, xlow, xhigh = h["xbins"]

                ybins = np.array([])

                for bs in h["ybins"]:
                    low, high, step = bs
                    ybins = np.append(ybins, np.arange(low, high+step, step))

                logging.debug(f"Creating 2D histogram {key} with variable y bins:\n{ybins} and {nbinsx} regular x bins from {xlow} to {xhigh}")

                hists.append(df_ref.Histo2D((hist_name, h["title"], nbinsx, xlow, xhigh, len(ybins)-1, ybins), h["xvar"], h["yvar"], weight))
            elif _is_variable_bins(h["xbins"]) and _is_variable_bins(h["ybins"]):
                xbins = np.array([])
                for bs in h["xbins"]:
                    low, high, step = bs
                    xbins = np.append(xbins, np.arange(low, high+step, step))

                ybins = np.array([])
                for bs in h["ybins"]:
                    low, high, step = bs
                    ybins = np.append(ybins, np.arange(low, high+step, step))

                logging.debug(f"Creating 2D histogram {key} with variable x and y bins:\n{xbins}\n and\n{ybins}")
                hists.append(df_ref.Histo2D((hist_name, h["title"], len(xbins)-1, xbins, len(ybins)-1, ybins), h["xvar"], h["yvar"], weight))
            else:
                raise ValueError(f"Invalid binning format for histogram {key}. Exiting...")

    save_branches: Optional[list[str]] = cfg.get("save_branches", None)

    if save_branches is not None:
        tree_name = f"fluxTree_{det_loc_str}_{sample_name}"
        tree_log_str = f"Preparing Tree '{tree_name}' with {len(save_branches)} branches:\n\n"
        for branch in save_branches:
            tree_log_str += "* " + branch + "\n"

        opts = ROOT.RDF.RSnapshotOptions()  # type: ignore
        opts.fMode = "UPDATE"

        if len(save_branches) == 1 and save_branches[0] == "*":
            logging.info(tree_log_str)
            logging.info(
                f"Snapshotting to {out_fname}. Event loop will be executed now, this might take a while..."
            )
            df.Snapshot(tree_name, out_fname, list(cfg["aliases"].keys()) + list(cfg["definitions"].keys())) #, opts)
        elif len(save_branches) > 0:
            logging.info(tree_log_str)
            logging.info(
                f"Snapshotting to {out_fname}. Event loop will be executed now, this might take a while..."
            )
            df.Snapshot(tree_name, out_fname, save_branches, opts)  # type: ignore
    else:
        logging.info("Not saving any branches.")

    open_fn: Callable[[str], uproot.WritableDirectory] = uproot.update if Path(out_fname).exists() else uproot.recreate

    with open_fn(out_fname) as f:
        for h in hists:
            hname = h.GetName()
            h_obj = h.GetValue()
            logging.info(f"Writing histogram {hname} to {out_fname}:{det_loc_str}/...")
            logging.debug(f"{h}, {h_obj}")
            f[det_loc_str + "/" + hname] = h_obj

        logging.info(f"Writing {pot} POT to {out_fname}:{sample_name}/hpot")

        hpot = ROOT.TH1D("hpot", "POT", 1, 0, 1)  # type: ignore
        hpot.SetBinContent(1, pot)

        f[f"{det_loc_str}/{sample_name}/hpot"] = hpot


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

    with open(args.config, "rb") as f:
        cfg = tomllib.load(f)

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

    for sample_name, files in cfg["file_sets"].items():
        parent_dir, pattern = files.rsplit("/", 1)

        file_list = [file for file in Path(parent_dir).glob(pattern) if file.is_file() and file.stat().st_size > 0]
        file_list = file_list[:cfg.get("max_files", None)]
        n_files = len(file_list)

        logging.info(f"Found {n_files} files in {parent_dir}")

        if n_files == 0:
            logging.error(f"Skipping {sample_name} as no files were found...")
            continue

        vec_files = ROOT.std.vector(ROOT.std.string)()
        for file in file_list:
            vec_files.push_back(str(file))
        if isinstance(cfg["location"], Iterable):
            det_locs = [[float(x) for x in loc] for loc in cfg["location"]]
        else:
            det_locs  = [[float(x) for x in cfg["location"]]]

        for det_loc in det_locs:
            if len(det_loc) != 3:
                logging.error(
                    f"Invalid location format: {det_loc}. Please provide '[x, y, z]'. Exiting..."
                )
                sys.exit(1)
            logging.info(f"Going to calculate weights for location: {det_loc}")

            logging.debug(f"Location: {det_loc}")

            # det_locs = np.array([
            #         np.arange(0, 1000, step=100),
            #         np.zeros(10),
            #         np.full(10, 10000)
            #         ])
            run_analysis(vec_files, str(out_fname), sample_name, det_loc, cfg, debug=args.debug)


if __name__ == "__main__":
    start = time.perf_counter()
    try:
        main()
    except KeyboardInterrupt:
        logging.info("Received keyboard interrupt. Exiting...")
        sys.exit(0)
    end = time.perf_counter()
    print(f"\nFinished in {datetime.timedelta(seconds=end - start)}")
