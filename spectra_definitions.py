import ROOT  # type: ignore


def apply_defs(df: ROOT.RDataFrame, det_loc: list[float]) -> ROOT.RDataFrame:
    _definitions = (
        df.Alias("nu_pdg", "dk2nu.decay.ntype")
        .Alias("parent_pdg", "dk2nu.decay.ptype")
        .Alias("vx", "dk2nu.decay.vx")
        .Alias("vy", "dk2nu.decay.vy")
        .Alias("vz", "dk2nu.decay.vz")
        .Alias("pdpx", "dk2nu.decay.pdpx")
        .Alias("pdpy", "dk2nu.decay.pdpy")
        .Alias("pdpz", "dk2nu.decay.pdpz")
        .Define("det_loc", f"ROOT::RVecD{{ {det_loc[0]}, {det_loc[1]}, {det_loc[2]} }}")
        .Define(
            "decay_vertex",
            "ROOT::RVecD{vx, vy, vz}",
        )
        .Define("rr", "det_loc - decay_vertex")
        .Define(
            "parent_momentum",
            "ROOT::RVecD{pdpx, pdpy, pdpz}",
        )
        .Define("sangdet", "Numba::calc_solid_angle(rr)")
        .Define("costh", "Numba::calc_costheta_par(parent_momentum, rr)")
        .Define("parent_mass", "Numba::pdg_to_mass(parent_pdg)")
        .Define("parent_energy", "Numba::calc_energy(parent_momentum, parent_mass)")
        .Define("pgamma", "Numba::calc_gamma(parent_energy, parent_mass)")
        .Define("emrat", "Numba::calc_energy_in_beam(pgamma, costh)")
        .Define("nu_energy", "dk2nu.decay.necm * emrat")
        .Define(
            "weight",
            "calc_weight(dk2nu.decay, sangdet, emrat, parent_energy, nu_energy, pgamma, rr)",
        )
        .Define("theta_p", "Numba::theta_p(parent_momentum, det_loc)")
        .Define("par_codes", "Numba::parent_to_code(dk2nu.ancestor.pdg)")
        .Define("target_codes", "Numba::target_to_code(dk2nu.ancestor.nucleus)")
    )
    return _definitions


_nu_pdg = {"numu": 14, "nue": 12, "numubar": -14, "nuebar": -12}


def build_histograms(horn: str, nu: str):
    energy_binning = [200, 0, 20]
    angle_binning = [180, -90, 90]

    nu_filter = f"dk2nu.decay.ntype == {_nu_pdg[nu]}"

    histograms = {
        f"hnu_E_{horn}_{nu}": lambda df: df.Filter(nu_filter).Histo1D(
            (f"hnu_E_{horn}_{nu}", "Nu Energy", *energy_binning), "nu_energy", "weight"
        ),
        f"hnu_theta_{horn}_{nu}": lambda df: df.Filter(nu_filter).Histo1D(
            (f"hnu_theta_{horn}_{nu}", "Parent Angle", *angle_binning),
            "theta_p",
            "weight",
        ),
        f"hnu_thetaE_{horn}_{nu}": lambda df: df.Filter(nu_filter).Histo2D(
            (
                f"hnu_thetaE_{horn}_{nu}",
                "Enu vs. Theta",
                *angle_binning,
                *energy_binning,
            ),
            "theta_p",
            "nu_energy",
            "weight",
        ),
        f"hancestors_{horn}_{nu}": lambda df: df.Filter(
            f"nu_energy > 0.400 && {nu_filter}"
        ).Histo2D(
            (f"hancestors_{horn}_{nu}", "Ancestor PDG", 5, 0, 5, 10, 0, 10),
            "target_codes",
            "par_codes",
            "weight",
        ),
    }
    return histograms
