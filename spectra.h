#ifndef SPECTRA_H
#define SPECTRA_H

#include "pdgid.h"

#include "dk2nu.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <string_view>
#include <unordered_map>
#include <utility>

using pdg::PDGID;

class Spectra {
public:
  Spectra(std::string_view id, const int nu_pdg);

  void FillSpectra(const bsim::Dk2Nu& dk2nu,
                   const double wght,
                   const double nu_energy,
                   const double theta_par);
  void WriteHistograms();

private:
  double nNeutrinos = 0.;
  double nInteractions = 0.;

  TH1D hnu_E;
  TH1D hnu_E_pipm;
  TH1D hnu_E_kpm;
  TH1D hnu_E_k0l;
  TH1D hnu_E_mu;

  TH1D hnu_theta;
  TH1D hnu_theta_pipm;
  TH1D hnu_theta_kpm;
  TH1D hnu_theta_k0l;
  TH1D hnu_theta_mu;

  TH1I hnInteractions;
  TH1I hnNeutrinos;
  TH2F hancestorInteractions;

  TH2D hnu_thetaE;
  TH2D hnu_thetaE_pipm;
  TH2D hnu_thetaE_kpm;
  TH2D hnu_thetaE_k0l;
  TH2D hnu_thetaE_mu;

  TH2D hnu_thetaY;
  TH2D hnu_thetaY_pipm;
  TH2D hnu_thetaY_kpm;
  TH2D hnu_thetaY_k0l;
  TH2D hnu_thetaY_mu;

  TH2D hnu_thetaZ;
  TH2D hnu_thetaZ_pipm;
  TH2D hnu_thetaZ_kpm;
  TH2D hnu_thetaZ_k0l;
  TH2D hnu_thetaZ_mu;

  TH2D hnu_zE;
  TH2D hnu_zE_pipm;
  TH2D hnu_zE_kpm;
  TH2D hnu_zE_k0l;
  TH2D hnu_zE_mu;

  TH2D hnu_yE;
  TH2D hnu_yE_pipm;
  TH2D hnu_yE_kpm;
  TH2D hnu_yE_k0l;
  TH2D hnu_yE_mu;

  TH2D hnu_xy;
  TH2D hnu_zx;
  TH2D hnu_zy;

  TH2D hnu_xy_pipm;
  TH2D hnu_zx_pipm;
  TH2D hnu_zy_pipm;

  TH2D hnu_xy_kpm;
  TH2D hnu_zx_kpm;
  TH2D hnu_zy_kpm;

  TH2D hnu_xy_k0l;
  TH2D hnu_zx_k0l;
  TH2D hnu_zy_k0l;

  TH2D hnu_xy_mu;
  TH2D hnu_zx_mu;
  TH2D hnu_zy_mu;

  TH2D hpar_xy;
  TH2D hpar_zx;
  TH2D hpar_zy;

  TH3D hnu_xyz;
  TH3D hpar_xyz;

  std::unordered_map<PDGID, float> pdg2Index;
  std::vector<std::pair<PDGID, PDGID>> IdentifyPrecursors(const bsim::Dk2Nu& dk2nu);
};

#endif
