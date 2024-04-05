#include "spectra.h"
#include "helper.h"
#include "vecops.h"

#include "TChain.h"
#include "TFile.h"
#include "TStopwatch.h"

using pdg::PDGID;
using helper::getFileList;

double GetWeight(const bsim::Dk2Nu& dk2nu,
                 const std::array<double, 3>& detCoords,
                 double& nu_energy)
{
  static constexpr double rdet = 100.0; //in cm

  const auto parent_ptype = static_cast<PDGID>(dk2nu.decay.ptype);

  if (pdg::pdgid2Mass.find(parent_ptype) == pdg::pdgid2Mass.end()) {
    std::cout << "GetWeight - Wrong parent type!! " << parent_ptype << '\n';
    return -9999.;
  }

  const double parent_mass = pdg::pdgid2Mass.at(parent_ptype);

  const double pdPx = dk2nu.decay.pdpx;
  const double pdPy = dk2nu.decay.pdpy;
  const double pdPz = dk2nu.decay.pdpz;

  const double parent_energy =
    std::sqrt((pdPx * pdPx) + (pdPy * pdPy) + (pdPz * pdPz) + (parent_mass * parent_mass));
  const double gamma = parent_energy / parent_mass;
  const double gamma_sqr = gamma * gamma;
  const double beta = std::sqrt((gamma_sqr - 1.) / gamma_sqr);

  const double parentP = std::sqrt((pdPx * pdPx) + (pdPy * pdPy) + (pdPz * pdPz));

  const double rrx = detCoords[0] - dk2nu.decay.vx;
  const double rry = detCoords[1] - dk2nu.decay.vy;
  const double rrz = detCoords[2] - dk2nu.decay.vz;
  const double rr = std::sqrt(rrx * rrx + rry * rry + rrz * rrz);

  double costh_parent = ((pdPx * rrx) + (pdPy * rry) + (pdPz * rrz)) / (parentP * rr);

  if (costh_parent > 1.) { costh_parent = 1.; }
  else if (costh_parent < -1.) {
    costh_parent = -1.;
  }
  const double emrat = 1. / (gamma * (1. - beta * costh_parent));
  const double angle_detector = (rdet * rdet) / (rr * rr) / 4.;

  double wght = (angle_detector * emrat * emrat * dk2nu.decay.nimpwt) / 3.1416;

  const double enuzr = dk2nu.decay.necm;

  nu_energy = emrat * enuzr;

  //done for all except polarized muon
  // in which case need to modify weight
  if (std::abs(parent_ptype) == PDGID::MUM) {
    //boost new neutrino to mu decay cm
    double beta[3];
    double p_nu[3]; //nu momentum
    beta[0] = pdPx / parent_energy;
    beta[1] = pdPy / parent_energy;
    beta[2] = pdPz / parent_energy;

    p_nu[0] = rrx * nu_energy / rr;
    p_nu[1] = rry * nu_energy / rr;
    p_nu[2] = rrz * nu_energy / rr;

    double partial = gamma * (beta[0] * p_nu[0] + beta[1] * p_nu[1] + beta[2] * p_nu[2]);
    partial = nu_energy - partial / (gamma + 1.);
    double p_dcm_nu[4];
    for (int i = 0; i < 3; i++)
      p_dcm_nu[i] = p_nu[i] - beta[i] * gamma * partial;
    p_dcm_nu[3] = 0.;
    for (int i = 0; i < 3; i++)
      p_dcm_nu[3] += p_dcm_nu[i] * p_dcm_nu[i];
    p_dcm_nu[3] = std::sqrt(p_dcm_nu[3]);

    //boost parent of mu to mu production cm
    // gamma = m_ppenergy / m_parent_mass;
    beta[0] = dk2nu.decay.ppdxdz * dk2nu.decay.pppz / dk2nu.decay.ppenergy;
    beta[1] = dk2nu.decay.ppdydz * dk2nu.decay.pppz / dk2nu.decay.ppenergy;
    beta[2] = dk2nu.decay.pppz / dk2nu.decay.ppenergy;
    partial = gamma * (beta[0] * dk2nu.decay.muparpx + beta[1] * dk2nu.decay.muparpy +
                       beta[2] * dk2nu.decay.muparpz);
    partial = dk2nu.decay.mupare - partial / (gamma + 1.);
    double p_pcm_mp[4];
    p_pcm_mp[0] = dk2nu.decay.muparpx - beta[0] * gamma * partial;
    p_pcm_mp[1] = dk2nu.decay.muparpy - beta[1] * gamma * partial;
    p_pcm_mp[2] = dk2nu.decay.muparpz - beta[2] * gamma * partial;
    p_pcm_mp[3] = 0.;
    for (int i = 0; i < 3; i++)
      p_pcm_mp[3] += p_pcm_mp[i] * p_pcm_mp[i];
    p_pcm_mp[3] = std::sqrt(p_pcm_mp[3]);

    double wt_ratio = 1.;
    //have to check p_pcm_mp
    //it can be 0 if mupar..=0. (I guess muons created in target??)
    if (p_pcm_mp[3] != 0.) {
      //calc new decay angle w.r.t. (anti)spin direction
      double costh =
        (p_dcm_nu[0] * p_pcm_mp[0] + p_dcm_nu[1] * p_pcm_mp[1] + p_dcm_nu[2] * p_pcm_mp[2]) /
        (p_dcm_nu[3] * p_pcm_mp[3]);

      if (costh > 1.)
        costh = 1.;
      else if (costh < -1.)
        costh = -1.;

      //calc relative weight due to angle difference
      auto const nu_type = dk2nu.decay.ntype;
      if (nu_type == PDGID::NUE || nu_type == PDGID::NUEBAR) { wt_ratio = 1. - costh; }
      else if (nu_type == PDGID::NUMU || nu_type == PDGID::NUMUBAR) {
        const double mumass = pdg::pdgid2Mass.at(PDGID::MUM);
        double xnu = 2. * enuzr / mumass;
        wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
      }
      else if (nu_type == PDGID::NUTAU || nu_type == PDGID::NUTAUBAR) {
        const double taumass = pdg::pdgid2Mass.at(PDGID::TAUM);
        double xnu = 2. * enuzr / taumass;
        wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
        std::cout << "calculating weight for tau neutrino; this may not be correct" << std::endl;
      }
      else {
        std::cout << "eventRates:: Bad neutrino type = " << nu_type << std::endl;
      }
    }
    wght *= wt_ratio;
  }

  return wght;
}

// static const std::array<double, 3> ICARUSCoords = {450.37, 7991.98, 79512.66};
static constexpr std::array<double, 3> TWOXTWOCoords = {0., 0., 103648.837};
static constexpr std::array<double, 3> NuMIAxis = {0., 0., 1.};

void runEventLoop(TChain& chain,
                  const bsim::Dk2Nu* dk2nu,
                  std::unordered_map<int, Spectra*>& spectra,
                  const unsigned int max_evts = 0)
{
  chain.SetBranchAddress("dk2nu", &dk2nu);

  auto const& coords = TWOXTWOCoords;

  static const auto n = vecops::cross(coords, NuMIAxis);
  static constexpr double magN = 1.;

  const unsigned int n_entries = chain.GetEntries();

  double pct = 0.;

  static constexpr double RAD2DEG = 180. / 3.1416;

  for (std::size_t entry = 0; entry < n_entries; ++entry) {
    if ((max_evts > 0) && (entry > max_evts)) { break; }
    pct = 100. * entry / n_entries;

    chain.GetEntry(entry);
    auto& spec = spectra[dk2nu->decay.ntype];
    double nu_energy = -9999.; // GeV
    const std::array<double, 3> parentDecayP = {
      dk2nu->decay.pdpx, dk2nu->decay.pdpy, dk2nu->decay.pdpz};
    const auto parentDecayPProj = vecops::project(parentDecayP, n);
    const auto magPProj = vecops::mag(parentDecayPProj);

    const double costheta_par = vecops::dot(parentDecayPProj, NuMIAxis) / (magPProj * magN);

    double theta_par = RAD2DEG * std::acos(costheta_par);

    // check if parentDecayPProj is between NuMIAxis and ICARUSCoords on the plane formed by NuMIAxis and ICARUSCoords
    // and negate theta if it is not.
    const auto crossProd = vecops::cross(parentDecayPProj, NuMIAxis);
    const double dotProd = vecops::dot(crossProd, n);

    if (dotProd < 0.) { theta_par *= -1.; }

    const double wght = GetWeight(*dk2nu, coords, nu_energy);

    spec->FillSpectra(*dk2nu, wght, nu_energy, theta_par);

    if ((entry == 0) || (entry % 10000 == 0)) {
      std::cout << "FILLING ENTRY " << entry << " / " << n_entries << " (" << pct << "%)"
                << " WITH WEIGHT = " << wght << '\n';
    }

  } // Event Loop
} // runEventLoop

void processFileList(const char* files, std::string_view fileIndex)
{
  auto const fileList = getFileList(files);

  auto pChain = std::make_unique<TChain>("dk2nuTree");

  for (auto const& f : fileList) {
    pChain->Add(f.c_str());
  }

  const bsim::Dk2Nu dk2nu;
  Spectra spec_fhc_numu(fileIndex, 14);
  Spectra spec_fhc_numubar(fileIndex, -14);
  Spectra spec_fhc_nue(fileIndex, 12);
  Spectra spec_fhc_nuebar(fileIndex, -12);

  std::unordered_map<int, Spectra*> spectra = {
    {14, &spec_fhc_numu}, {-14, &spec_fhc_numubar}, {12, &spec_fhc_nue}, {-12, &spec_fhc_nuebar}};

  runEventLoop(*pChain, &dk2nu, spectra);

  constexpr double pot_per_file = 500000.;
  const double total_pot = pot_per_file * static_cast<double>(fileList.size());

  const std::string potName = "pot_" + std::string(fileIndex);
  TH1D hpot(potName.c_str(), potName.c_str(), 1, 0, 1);
  hpot.SetBinContent(1, total_pot);

  hpot.Write();

  for (auto const& spec : spectra) {
    spec.second->WriteHistograms();
  }
}

void testPrint(std::string_view x, std::string_view y)
{
  std::cout << "Processing file list " << x << " with name " << y << '\n';
}

int main(int argc, char** argv)
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <output file> <file list 1> <name 1> <file list 2> <name 2> ... <file list n> "
                 "<name n>\n";
    return 1;
  }

  TStopwatch sw;
  sw.Start();

  const char* outName = argv[1];
  TFile outputFile(outName, "RECREATE");

  for (int i = 2; i < argc - 1; i += 2) {
    std::cout << "Processing file list " << argv[i] << " with name " << argv[i + 1] << '\n';
    processFileList(argv[i], std::string(argv[i + 1]));
  }

  outputFile.Close();

  std::cout << "Done in " << sw.RealTime() << " seconds." << std::endl;
  return 0;
}
