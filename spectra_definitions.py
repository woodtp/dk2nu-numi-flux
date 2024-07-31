import ROOT  # type: ignore


def apply_defs(df: ROOT.RDataFrame, det_loc: list[float]) -> ROOT.RDataFrame:  # type: ignore
    if len(det_loc) != 3:
        raise ValueError("Detector location must be a list of 3 floats!")

    _definitions = (
        df.Alias("nu_pdg", "dk2nu.decay.ntype")
        .Alias("parent_pdg", "dk2nu.decay.ptype")
        .Alias("vx", "dk2nu.decay.vx")
        .Alias("vy", "dk2nu.decay.vy")
        .Alias("vz", "dk2nu.decay.vz")
        .Alias("pdpx", "dk2nu.decay.pdpx")
        .Alias("pdpy", "dk2nu.decay.pdpy")
        .Alias("pdpz", "dk2nu.decay.pdpz")
        # .Define("is_new_g4", 'dk2nu.ancestor[0].proc == "BeamParticle"')
        .Define("ancestor_vol", "get_volumes(dk2nu.ancestor)")
        .Define(
            "det_loc",
            f"ROOT::RVec<double>{{ {det_loc[0]}, {det_loc[1]}, {det_loc[2]} }}",
        )
        .Define(
            "decay_vertex",
            "ROOT::RVec<double>{vx, vy, vz}",
        )
        .Define("rr", "det_loc - decay_vertex")
        .Define("parent_momentum", "sqrt(pdpx*pdpx + pdpy*pdpy + pdpz*pdpz)")
        .Define(
            "v_parent_momentum",
            "ROOT::RVec<double>{pdpx, pdpy, pdpz}",
        )
        .Define("sangdet", "Numba::calc_solid_angle(rr)")                               # Solid angle as fraction of 4pi
        .Define("costh", "Numba::calc_costheta_par(v_parent_momentum, rr)")             # Cosine of the angle between the parent particle and the detector component
        .Define("parent_mass", "Numba::pdg_to_mass(parent_pdg)")
        .Define("parent_energy", "Numba::calc_energy(v_parent_momentum, parent_mass)")  # Calculate the energy of the parent particle at the decay vertex
        .Define("pgamma", "Numba::calc_gamma(parent_energy, parent_mass)")              # Calculate the gamma factor of the parent particle
        .Define("emrat", "Numba::calc_energy_in_beam(pgamma, costh)")                   # Weighted neutrino energy in the lab frame. Small angle approximation.
        .Define("nu_energy", "dk2nu.decay.necm * emrat")                                # Lorentz boost the neutrino energy into the lab frame
        .Define(
            "weight",
            "calc_weight(dk2nu.decay, sangdet, emrat, parent_energy, nu_energy, pgamma, rr)",
        )
        .Define("theta_p", "Numba::theta_p(v_parent_momentum, det_loc)")
        .Define("par_codes", "Numba::parent_to_code(dk2nu.ancestor.pdg)")
        .Define("target_codes", "Numba::target_to_code(dk2nu.ancestor.nucleus)")
        .Define("ancestor_parent_pdg", "Numba::ancestor_parent_pdg(dk2nu.ancestor.pdg)")
        # .Define("ancestor_parent_mom", "Numba::calc_magnitudes(dk2nu.ancestor.pprodpx, dk2nu.ancestor.pprodpy, dk2nu.ancestor.pprodpz)")
        .Define("ancestor_parent_mom", "get_incident_momenta(dk2nu.ancestor)")
        .Define("ancestor_produced_mom", "get_produced_momenta(dk2nu.ancestor)")
        .Define("ancestor_produced_theta", "calc_theta(dk2nu.ancestor)")
        .Define("ancestor_mass", "Numba::ancestor_pdg2mass(dk2nu.ancestor.pdg)")
        .Define("ancestor_pT", "calc_pT(dk2nu.ancestor)")
        .Define("ancestor_xF", "calc_xF(dk2nu.ancestor, ancestor_mass)")
    )
    return _definitions
