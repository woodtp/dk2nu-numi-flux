#include "spectra.h"
#include <iostream>

Spectra::Spectra(std::string_view id, const int nu_pdg = -1)
{
  pdg2Index = {{PDGID::PROTON, 0.5},
               {PDGID::NEUTRON, 1.5},
               {PDGID::PIP, 2.5},
               {PDGID::PIM, 3.5},
               {PDGID::KP, 4.5},
               {PDGID::KM, 5.5},
               {PDGID::K0L, 6.5},
               {PDGID::MUM, 7.5},
               {PDGID::MUP, 8.5},
               {PDGID::CARBON, 0.5},
               {PDGID::ALUMINUM, 1.5},
               {PDGID::IRON, 2.5},
               {PDGID::NULLTARGET, 4.5},
               {PDGID::ZERO, 4.5}};

  static constexpr int nbinsx = 80;
  static constexpr double xlow = -200.;
  static constexpr double xup = 200.;

  static constexpr int nbinsy = 80;
  static constexpr double ylow = -200.;
  static constexpr double yup = 200.;

  static constexpr int nbinsz = 260;
  static constexpr double zlow = -200.;
  static constexpr double zup = 5000.;

  static constexpr int nbinsE = 200;
  static constexpr double elow = 0.;
  static constexpr double eup = 20.;

  static constexpr int nbinsTh = 180;
  static constexpr double thlow = -90.;
  static constexpr double thup = 90.;

  std::string nu_suff = std::abs(nu_pdg) == 12 ? "nue" : std::abs(nu_pdg) == 14 ? "numu" : "";
  if (nu_pdg < 0 && nu_pdg != -1) { nu_suff += "bar"; }

  const std::string suffix = std::string(id) + "_" + nu_suff;

  const std::string ints_label = "hancestorInteractions_" + suffix;

  const std::string elabel = "hnu_E_" + suffix;
  const std::string elabel_pipm = "hnu_E_pipm_" + suffix;
  const std::string elabel_kpm = "hnu_E_kpm_" + suffix;
  const std::string elabel_k0l = "hnu_E_k0l_" + suffix;
  const std::string elabel_mu = "hnu_E_mu_" + suffix;

  const std::string thetalabel = "hnu_theta_" + suffix;
  const std::string thetalabel_pipm = "hnu_theta_pipm_" + suffix;
  const std::string thetalabel_kpm = "hnu_theta_kpm_" + suffix;
  const std::string thetalabel_k0l = "hnu_theta_k0l_" + suffix;
  const std::string thetalabel_mu = "hnu_theta_mu_" + suffix;

  const std::string thetaElabel = "hnu_thetaE_" + suffix;
  const std::string thetaElabel_pipm = "hnu_thetaE_pipm_" + suffix;
  const std::string thetaElabel_kpm = "hnu_thetaE_kpm_" + suffix;
  const std::string thetaElabel_k0l = "hnu_thetaE_k0l_" + suffix;
  const std::string thetaElabel_mu = "hnu_thetaE_mu_" + suffix;

  const std::string thetaYlabel = "hnu_thetaY_" + suffix;
  const std::string thetaYlabel_pipm = "hnu_thetaY_pipm_" + suffix;
  const std::string thetaYlabel_kpm = "hnu_thetaY_kpm_" + suffix;
  const std::string thetaYlabel_k0l = "hnu_thetaY_k0l_" + suffix;
  const std::string thetaYlabel_mu = "hnu_thetaY_mu_" + suffix;

  const std::string thetaZlabel = "hnu_thetaZ_" + suffix;
  const std::string thetaZlabel_pipm = "hnu_thetaZ_pipm_" + suffix;
  const std::string thetaZlabel_kpm = "hnu_thetaZ_kpm_" + suffix;
  const std::string thetaZlabel_k0l = "hnu_thetaZ_k0l_" + suffix;
  const std::string thetaZlabel_mu = "hnu_thetaZ_mu_" + suffix;

  const std::string zElabel = "hnu_zE_" + suffix;
  const std::string zElabel_pipm = "hnu_zE_pipm_" + suffix;
  const std::string zElabel_kpm = "hnu_zE_kpm_" + suffix;
  const std::string zElabel_k0l = "hnu_zE_k0l_" + suffix;
  const std::string zElabel_mu = "hnu_zE_mu_" + suffix;

  const std::string yElabel = "hnu_yE_" + suffix;
  const std::string yElabel_pipm = "hnu_yE_pipm_" + suffix;
  const std::string yElabel_kpm = "hnu_yE_kpm_" + suffix;
  const std::string yElabel_k0l = "hnu_yE_k0l_" + suffix;
  const std::string yElabel_mu = "hnu_yE_mu_" + suffix;

  const std::string zylabel = "hnu_zy_" + suffix;
  const std::string zylabel_pipm = "hnu_zy_pipm_" + suffix;
  const std::string zylabel_kpm = "hnu_zy_kpm_" + suffix;
  const std::string zylabel_k0l = "hnu_zy_k0l_" + suffix;
  const std::string zylabel_mu = "hnu_zy_mu_" + suffix;

  const std::string zxlabel = "hnu_zx_" + suffix;
  const std::string zxlabel_pipm = "hnu_zx_pipm_" + suffix;
  const std::string zxlabel_kpm = "hnu_zx_kpm_" + suffix;
  const std::string zxlabel_k0l = "hnu_zx_k0l_" + suffix;
  const std::string zxlabel_mu = "hnu_zx_mu_" + suffix;

  const std::string xylabel = "hnu_xy_" + suffix;
  const std::string xylabel_pipm = "hnu_xy_pipm_" + suffix;
  const std::string xylabel_kpm = "hnu_xy_kpm_" + suffix;
  const std::string xylabel_k0l = "hnu_xy_k0l_" + suffix;
  const std::string xylabel_mu = "hnu_xy_mu_" + suffix;

  const std::string zylabel_par = "hpar_zy_" + suffix;
  const std::string zxlabel_par = "hpar_zx_" + suffix;
  const std::string xylabel_par = "hpar_xy_" + suffix;

  const std::string xyzlabel = "hnu_xyz_" + suffix;
  const std::string xyzlabel_par = "hpar_xyz_" + suffix;

  const std::string nints_label = "hnInteractions_" + suffix;
  const std::string nnus_label = "hnNus_" + suffix;

  hnInteractions = TH1I(nints_label.c_str(), nints_label.c_str(), 1, 0., 1.);
  hnNeutrinos = TH1I(nnus_label.c_str(), nnus_label.c_str(), 1, 0., 1.);
  hancestorInteractions = TH2F(ints_label.c_str(), ints_label.c_str(), 5, 0., 5., 10, 0., 10.);

  hnu_E = TH1D(elabel.c_str(), elabel.c_str(), nbinsE, elow, eup);
  hnu_E_pipm = TH1D(elabel_pipm.c_str(), elabel_pipm.c_str(), nbinsE, elow, eup);
  hnu_E_kpm = TH1D(elabel_kpm.c_str(), elabel_kpm.c_str(), nbinsE, elow, eup);
  hnu_E_k0l = TH1D(elabel_k0l.c_str(), elabel_k0l.c_str(), nbinsE, elow, eup);
  hnu_E_mu = TH1D(elabel_mu.c_str(), elabel_mu.c_str(), nbinsE, elow, eup);

  hnu_theta = TH1D(thetalabel.c_str(), thetalabel.c_str(), nbinsTh, thlow, thup);
  hnu_theta_pipm = TH1D(thetalabel_pipm.c_str(), thetalabel_pipm.c_str(), nbinsTh, thlow, thup);
  hnu_theta_kpm = TH1D(thetalabel_kpm.c_str(), thetalabel_kpm.c_str(), nbinsTh, thlow, thup);
  hnu_theta_k0l = TH1D(thetalabel_k0l.c_str(), thetalabel_k0l.c_str(), nbinsTh, thlow, thup);
  hnu_theta_mu = TH1D(thetalabel_mu.c_str(), thetalabel_mu.c_str(), nbinsTh, thlow, thup);

  hnu_thetaE =
    TH2D(thetaElabel.c_str(), thetaElabel.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
  hnu_thetaE_pipm = TH2D(
    thetaElabel_pipm.c_str(), thetaElabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
  hnu_thetaE_kpm =
    TH2D(thetaElabel_kpm.c_str(), thetaElabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
  hnu_thetaE_k0l =
    TH2D(thetaElabel_k0l.c_str(), thetaElabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
  hnu_thetaE_mu =
    TH2D(thetaElabel_mu.c_str(), thetaElabel_mu.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);

  hnu_thetaY =
    TH2D(thetaYlabel.c_str(), thetaYlabel.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
  hnu_thetaY_pipm = TH2D(
    thetaYlabel_pipm.c_str(), thetaYlabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
  hnu_thetaY_kpm =
    TH2D(thetaYlabel_kpm.c_str(), thetaYlabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
  hnu_thetaY_k0l =
    TH2D(thetaYlabel_k0l.c_str(), thetaYlabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
  hnu_thetaY_mu =
    TH2D(thetaYlabel_mu.c_str(), thetaYlabel_mu.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);

  hnu_thetaZ =
    TH2D(thetaZlabel.c_str(), thetaZlabel.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
  hnu_thetaZ_pipm = TH2D(
    thetaZlabel_pipm.c_str(), thetaZlabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
  hnu_thetaZ_kpm =
    TH2D(thetaZlabel_kpm.c_str(), thetaZlabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
  hnu_thetaZ_k0l =
    TH2D(thetaZlabel_k0l.c_str(), thetaZlabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
  hnu_thetaZ_mu =
    TH2D(thetaZlabel_mu.c_str(), thetaZlabel_mu.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);

  hnu_zE = TH2D(zElabel.c_str(), zElabel.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
  hnu_zE_pipm =
    TH2D(zElabel_pipm.c_str(), zElabel_pipm.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
  hnu_zE_kpm = TH2D(zElabel_kpm.c_str(), zElabel_kpm.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
  hnu_zE_k0l = TH2D(zElabel_k0l.c_str(), zElabel_k0l.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
  hnu_zE_mu = TH2D(zElabel_mu.c_str(), zElabel_mu.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);

  hnu_yE = TH2D(yElabel.c_str(), yElabel.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
  hnu_yE_pipm =
    TH2D(yElabel_pipm.c_str(), yElabel_pipm.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
  hnu_yE_kpm = TH2D(yElabel_kpm.c_str(), yElabel_kpm.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
  hnu_yE_k0l = TH2D(yElabel_k0l.c_str(), yElabel_k0l.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
  hnu_yE_mu = TH2D(yElabel_mu.c_str(), yElabel_mu.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);

  hnu_zy = TH2D(zylabel.c_str(), zylabel.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
  hnu_zx = TH2D(zxlabel.c_str(), zxlabel.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
  hnu_xy = TH2D(xylabel.c_str(), xylabel.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

  hnu_zy_pipm =
    TH2D(zylabel_pipm.c_str(), zylabel_pipm.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
  hnu_zx_pipm =
    TH2D(zxlabel_pipm.c_str(), zxlabel_pipm.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
  hnu_xy_pipm =
    TH2D(xylabel_pipm.c_str(), xylabel_pipm.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

  hnu_zy_kpm = TH2D(zylabel_kpm.c_str(), zylabel_kpm.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
  hnu_zx_kpm = TH2D(zxlabel_kpm.c_str(), zxlabel_kpm.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
  hnu_xy_kpm = TH2D(xylabel_kpm.c_str(), xylabel_kpm.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

  hnu_zy_k0l = TH2D(zylabel_k0l.c_str(), zylabel_k0l.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
  hnu_zx_k0l = TH2D(zxlabel_k0l.c_str(), zxlabel_k0l.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
  hnu_xy_k0l = TH2D(xylabel_k0l.c_str(), xylabel_k0l.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

  hnu_zy_mu = TH2D(zylabel_mu.c_str(), zylabel_mu.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
  hnu_zx_mu = TH2D(zxlabel_mu.c_str(), zxlabel_mu.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
  hnu_xy_mu = TH2D(xylabel_mu.c_str(), xylabel_mu.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

  hpar_zy = TH2D(zylabel_par.c_str(), zylabel_par.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
  hpar_zx = TH2D(zxlabel_par.c_str(), zxlabel_par.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
  hpar_xy = TH2D(xylabel_par.c_str(), xylabel_par.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

  hnu_xyz = TH3D(
    xyzlabel.c_str(), xyzlabel.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
  hpar_xyz = TH3D(xyzlabel_par.c_str(),
                  xyzlabel_par.c_str(),
                  nbinsx,
                  xlow,
                  xup,
                  nbinsy,
                  ylow,
                  yup,
                  nbinsz,
                  zlow,
                  zup);
};

std::vector<std::pair<PDGID, PDGID>> Spectra::IdentifyPrecursors(const bsim::Dk2Nu& dk2nu)
{
  std::vector<std::pair<PDGID, PDGID>> precursors;
  precursors.reserve(dk2nu.ancestor.size());
  const bool isUpdated =
    dk2nu.ancestor[0].nucleus == 0; // account for the change in convention between versions
  for (std::size_t i = 0; i < dk2nu.ancestor.size() - 1; ++i) {
    auto const tgt_idx = isUpdated ? i + 1 : i;
    const auto parent_ptype = static_cast<PDGID>(dk2nu.ancestor[i].pdg);
    const auto target_type = static_cast<PDGID>(dk2nu.ancestor[tgt_idx].nucleus);
    precursors.emplace_back(std::make_pair(parent_ptype, target_type));
  }
  return precursors;
}

void Spectra::FillSpectra(const bsim::Dk2Nu& dk2nu,
                          const double wght,
                          const double nu_energy,
                          const double theta_par)
{

  const double nu_vx = dk2nu.decay.vx;
  const double nu_vy = dk2nu.decay.vy;
  const double nu_vz = dk2nu.decay.vz;

  const double par_vx = dk2nu.ppvx;
  const double par_vy = dk2nu.ppvy;
  const double par_vz = dk2nu.ppvz;

  auto const ptype = std::abs(dk2nu.decay.ptype);

  if (nu_energy > 0.4) {
    nNeutrinos += wght;

    auto const& precursors = IdentifyPrecursors(dk2nu);

    nInteractions += wght * precursors.size();

    // int i = 1;
    for (auto const& precursor : precursors) {
      auto const& target = precursor.second;
      auto const& proj = precursor.first;

      auto const projIdx = pdg2Index.find(proj) != pdg2Index.end() ? pdg2Index.at(proj) : 9.5;
      auto const targetIdx = pdg2Index.find(target) != pdg2Index.end() ? pdg2Index.at(target) : 3.5;
      // if (proj == PDGID::MUP || proj == PDGID::MUM) {
      //   std::cout << i << " / " << precursors.size() << " >> " << proj << " " << target << std::endl;
      //   std::cout << "    Filling " << projIdx << " " << targetIdx << std::endl;
      //   // continue;
      // }
      hancestorInteractions.Fill(targetIdx, projIdx, wght);
      // i++;
    }
  }

  hnu_E.Fill(nu_energy, wght);
  hnu_theta.Fill(theta_par, wght);
  hnu_thetaE.Fill(theta_par, nu_energy, wght);
  hnu_thetaY.Fill(theta_par, nu_vy, wght);
  hnu_thetaZ.Fill(theta_par, nu_vz, wght);
  hnu_yE.Fill(nu_vy, nu_energy, wght);
  hnu_zE.Fill(nu_vz, nu_energy, wght);

  hnu_zy.Fill(nu_vz, nu_vy, wght);
  hnu_zx.Fill(nu_vz, nu_vx, wght);
  hnu_xy.Fill(nu_vx, nu_vy, wght);
  hnu_xyz.Fill(nu_vx, nu_vy, nu_vz, wght);

  hpar_zy.Fill(par_vz, par_vy, wght);
  hpar_zx.Fill(par_vz, par_vx, wght);
  hpar_xy.Fill(par_vx, par_vy, wght);
  hpar_xyz.Fill(par_vx, par_vy, par_vz, wght);

  if (ptype == PDGID::PIP) {
    hnu_E_pipm.Fill(nu_energy, wght);
    hnu_zy_pipm.Fill(nu_vz, nu_vy, wght);
    hnu_zx_pipm.Fill(nu_vz, nu_vx, wght);
    hnu_xy_pipm.Fill(nu_vx, nu_vy, wght);
    hnu_theta_pipm.Fill(theta_par, wght);
    hnu_thetaE_pipm.Fill(theta_par, nu_energy, wght);
    hnu_thetaY_pipm.Fill(theta_par, nu_vy, wght);
    hnu_thetaZ_pipm.Fill(theta_par, nu_vz, wght);
    hnu_yE_pipm.Fill(nu_vy, nu_energy, wght);
    hnu_zE_pipm.Fill(nu_vz, nu_energy, wght);
  }
  else if (ptype == PDGID::KP) {
    hnu_E_kpm.Fill(nu_energy, wght);
    hnu_zy_kpm.Fill(nu_vz, nu_vy, wght);
    hnu_zx_kpm.Fill(nu_vz, nu_vx, wght);
    hnu_xy_kpm.Fill(nu_vx, nu_vy, wght);
    hnu_theta_kpm.Fill(theta_par, wght);
    hnu_thetaE_kpm.Fill(theta_par, nu_energy, wght);
    hnu_thetaY_kpm.Fill(theta_par, nu_vy, wght);
    hnu_thetaZ_kpm.Fill(theta_par, nu_vz, wght);
    hnu_yE_kpm.Fill(nu_vy, nu_energy, wght);
    hnu_zE_kpm.Fill(nu_vz, nu_energy, wght);
  }
  else if (ptype == PDGID::K0L) {
    hnu_E_k0l.Fill(nu_energy, wght);
    hnu_zy_k0l.Fill(nu_vz, nu_vy, wght);
    hnu_zx_k0l.Fill(nu_vz, nu_vx, wght);
    hnu_xy_k0l.Fill(nu_vx, nu_vy, wght);
    hnu_theta_k0l.Fill(theta_par, wght);
    hnu_thetaE_k0l.Fill(theta_par, nu_energy, wght);
    hnu_thetaY_k0l.Fill(theta_par, nu_vy, wght);
    hnu_thetaZ_k0l.Fill(theta_par, nu_vz, wght);
    hnu_yE_k0l.Fill(nu_vy, nu_energy, wght);
    hnu_zE_k0l.Fill(nu_vz, nu_energy, wght);
  }
  else if (ptype == PDGID::MUM) {
    hnu_E_mu.Fill(nu_energy, wght);
    hnu_zy_mu.Fill(nu_vz, nu_vy, wght);
    hnu_zx_mu.Fill(nu_vz, nu_vx, wght);
    hnu_xy_mu.Fill(nu_vx, nu_vy, wght);
    hnu_theta_mu.Fill(theta_par, wght);
    hnu_thetaE_mu.Fill(theta_par, nu_energy, wght);
    hnu_thetaY_mu.Fill(theta_par, nu_vy, wght);
    hnu_thetaZ_mu.Fill(theta_par, nu_vz, wght);
    hnu_yE_mu.Fill(nu_vy, nu_energy, wght);
    hnu_zE_mu.Fill(nu_vz, nu_energy, wght);
  }
}

void Spectra::WriteHistograms()
{
  std::cout << "TOTAL NUMBER NUS: " << nNeutrinos << std::endl;
  hnNeutrinos.SetBinContent(1, nNeutrinos);
  hnInteractions.SetBinContent(1, nInteractions);

  // Compute the average number of interactions per neutrino
  hancestorInteractions.Scale(1. / nNeutrinos);

  hnNeutrinos.Write();
  hnInteractions.Write();
  hancestorInteractions.Write();

  hnu_E.Write();
  hnu_E_pipm.Write();
  hnu_E_kpm.Write();
  hnu_E_k0l.Write();
  hnu_E_mu.Write();

  hnu_theta.Write();
  hnu_theta_pipm.Write();
  hnu_theta_kpm.Write();
  hnu_theta_k0l.Write();
  hnu_theta_mu.Write();

  hnu_thetaE.Write();
  hnu_thetaE_pipm.Write();
  hnu_thetaE_kpm.Write();
  hnu_thetaE_k0l.Write();
  hnu_thetaE_mu.Write();

  hnu_thetaY.Write();
  hnu_thetaY_pipm.Write();
  hnu_thetaY_kpm.Write();
  hnu_thetaY_k0l.Write();
  hnu_thetaY_mu.Write();

  hnu_thetaZ.Write();
  hnu_thetaZ_pipm.Write();
  hnu_thetaZ_kpm.Write();
  hnu_thetaZ_k0l.Write();
  hnu_thetaZ_mu.Write();

  hnu_zE.Write();
  hnu_zE_pipm.Write();
  hnu_zE_kpm.Write();
  hnu_zE_k0l.Write();
  hnu_zE_mu.Write();

  hnu_yE.Write();
  hnu_yE_pipm.Write();
  hnu_yE_kpm.Write();
  hnu_yE_k0l.Write();
  hnu_yE_mu.Write();

  hnu_xy.Write();
  hnu_zx.Write();
  hnu_zy.Write();

  hnu_xy_pipm.Write();
  hnu_zx_pipm.Write();
  hnu_zy_pipm.Write();

  hnu_xy_kpm.Write();
  hnu_zx_kpm.Write();
  hnu_zy_kpm.Write();

  hnu_xy_k0l.Write();
  hnu_zx_k0l.Write();
  hnu_zy_k0l.Write();

  hnu_xy_mu.Write();
  hnu_zx_mu.Write();
  hnu_zy_mu.Write();

  hpar_xy.Write();
  hpar_zx.Write();
  hpar_zy.Write();

  hnu_xyz.Write();
  hpar_xyz.Write();
}
