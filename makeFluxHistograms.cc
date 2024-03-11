#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

std::vector<std::string> getFileList(const char *file_name)
{
    std::vector<std::string> listOfFiles;

    std::ifstream fileList(file_name);

    std::string line;

    if (!fileList.is_open())
    {
        std::cout << "unable to open file " << file_name << "\n";
        return listOfFiles;
    }

    while (std::getline(fileList, line))
    {
        listOfFiles.push_back(line);
    }

    return listOfFiles;
}

double GetWeight(const bsim::Dk2Nu &dk2nu, const std::array<double, 3> &detCoords, double &nu_energy, double &theta_parent)
{
    const double rdet = 100.0;     //in cm
    const double pimass = 0.13957; //in GeV
    const double kmass = 0.49368;
    const double k0mass = 0.49767;
    const double mumass = 0.105658389;
    const double taumass = 1.77682;
    const int nue = 12;
    const int nuebar = -12;
    const int numu = 14;
    const int numubar = -14;
    const int nutau = 16;
    const int nutaubar = -16;
    const int muplus = -13;
    const int muminus = 13;

    double parent_mass = -9999.;
    const auto parent_ptype = dk2nu.decay.ptype;

    if (std::abs(parent_ptype) == 211)
        parent_mass = pimass;
    else if (std::abs(parent_ptype) == 321)
        parent_mass = kmass;
    else if (parent_ptype == 311 || parent_ptype == 130 || parent_ptype == 310)
        parent_mass = k0mass;
    else if (std::abs(parent_ptype) == 13)
        parent_mass = mumass;
    else
    {
        std::cout << "eventRates::GetWeight - Wrong parent type!! " << parent_ptype << " = " << parent_ptype << std::endl;
        return -9999.;
    }

    const double pdPx = dk2nu.decay.pdpx;
    const double pdPy = dk2nu.decay.pdpy;
    const double pdPz = dk2nu.decay.pdpz;

    const double parent_energy = std::sqrt((pdPx * pdPx) + (pdPy * pdPy) + (pdPz * pdPz) + (parent_mass * parent_mass));
    const double gamma = parent_energy / parent_mass;
    const double gamma_sqr = gamma * gamma;
    const double beta = std::sqrt((gamma_sqr - 1.) / gamma_sqr);

    const double parentP = std::sqrt((pdPx * pdPx) + (pdPy * pdPy) + (pdPz * pdPz));

    const double detx = detCoords[0];
    const double dety = detCoords[1];
    const double detz = detCoords[2];

    const double rrx = detx - dk2nu.decay.vx;
    const double rry = dety - dk2nu.decay.vy;
    const double rrz = detz - dk2nu.decay.vz;
    const double rr = std::sqrt(rrx * rrx + rry * rry + rrz * rrz);

    double costh_parent = ((pdPx * rrx) + (pdPy * rry) + (pdPz * rrz)) / (parentP * rr);

    if (costh_parent > 1.)
        costh_parent = 1.;
    else if (costh_parent < -1.)
        costh_parent = -1.;
    theta_parent = std::acos(costh_parent);
    const double emrat = 1. / (gamma * (1. - beta * costh_parent));
    const double angle_detector = (rdet * rdet) / (rr * rr) / 4.;

    double wght = (angle_detector * emrat * emrat * dk2nu.decay.nimpwt) / 3.1416;

    const double enuzr = dk2nu.decay.necm;

    nu_energy = emrat * enuzr;

    //done for all except polarized muon
    // in which case need to modify weight
    if (parent_ptype == muplus || parent_ptype == muminus)
    {
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
        partial = gamma * (beta[0] * dk2nu.decay.muparpx + beta[1] * dk2nu.decay.muparpy + beta[2] * dk2nu.decay.muparpz);
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
        if (p_pcm_mp[3] != 0.)
        {
            //calc new decay angle w.r.t. (anti)spin direction
            double costh = (p_dcm_nu[0] * p_pcm_mp[0] + p_dcm_nu[1] * p_pcm_mp[1] + p_dcm_nu[2] * p_pcm_mp[2]) / (p_dcm_nu[3] * p_pcm_mp[3]);

            if (costh > 1.)
                costh = 1.;
            else if (costh < -1.)
                costh = -1.;

            //calc relative weight due to angle difference
            auto const nu_type = dk2nu.decay.ntype;
            if (nu_type == nue || nu_type == nuebar)
            {
                wt_ratio = 1. - costh;
            }
            else if (nu_type == numu || nu_type == numubar)
            {
                double xnu = 2. * enuzr / mumass;
                wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
            }
            else if (nu_type == nutau || nu_type == nutaubar)
            {
                double xnu = 2. * enuzr / taumass;
                wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
                std::cout << "calculating weight for tau neutrino; this may not be correct" << std::endl;
            }
            else
            {
                std::cout << "eventRates:: Bad neutrino type = " << nu_type << std::endl;
            }
        }
        wght *= wt_ratio;
    }

    return wght;
}

class Spectra
{
public:
    Spectra(const std::string& hist_name)
    {
        const int nbinsx = 200;
        const double xlow = -200.;
        const double xup = 200.;

        const int nbinsy = 200;
        const double ylow = -200.;
        const double yup = 200.;

        const int nbinsz = 520;
        const double zlow = -200.;
        const double zup = 5000.;

        const int nbinsE = 1000;
        const double elow = 0.;
        const double eup = 20.;

        // in mrad
        const int nbinsTh = 1000;
        const double thlow = 0.;
        const double thup = 10.;

        const std::string elabel = "hnu_E_" + hist_name;
        const std::string elabel_pipm = "hnu_E_pipm_" + hist_name;
        const std::string elabel_kpm = "hnu_E_kpm_" + hist_name;
        const std::string elabel_k0l = "hnu_E_k0l_" + hist_name;

        const std::string thetalabel_pipm = "hnu_theta_pipm_" + hist_name;
        const std::string thetalabel_kpm = "hnu_theta_kpm_" + hist_name;
        const std::string thetalabel_k0l = "hnu_theta_k0l_" + hist_name;

        const std::string thetaElabel_pipm = "hnu_thetaE_pipm_" + hist_name;
        const std::string thetaElabel_kpm = "hnu_thetaE_kpm_" + hist_name;
        const std::string thetaElabel_k0l = "hnu_thetaE_k0l_" + hist_name;

        const std::string zElabel_pipm = "hnu_zE_pipm_" + hist_name;
        const std::string zElabel_kpm = "hnu_zE_kpm_" + hist_name;
        const std::string zElabel_k0l = "hnu_zE_k0l_" + hist_name;

        const std::string yElabel_pipm = "hnu_yE_pipm_" + hist_name;
        const std::string yElabel_kpm = "hnu_yE_kpm_" + hist_name;
        const std::string yElabel_k0l = "hnu_yE_k0l_" + hist_name;

        const std::string zylabel = "hnu_zy_" + hist_name;
        const std::string zylabel_pipm = "hnu_zy_pipm_" + hist_name;
        const std::string zylabel_kpm = "hnu_zy_kpm_" + hist_name;
        const std::string zylabel_k0l = "hnu_zy_k0l_" + hist_name;

        const std::string zxlabel = "hnu_zx_" + hist_name;
        const std::string zxlabel_pipm = "hnu_zx_pipm_" + hist_name;
        const std::string zxlabel_kpm = "hnu_zx_kpm_" + hist_name;
        const std::string zxlabel_k0l = "hnu_zx_k0l_" + hist_name;

        const std::string xylabel = "hnu_xy_" + hist_name;
        const std::string xylabel_pipm = "hnu_xy_pipm_" + hist_name;
        const std::string xylabel_kpm = "hnu_xy_kpm_" + hist_name;
        const std::string xylabel_k0l = "hnu_xy_k0l_" + hist_name;

        const std::string zylabel_par = "hpar_zy_" + hist_name;
        const std::string zxlabel_par = "hpar_zx_" + hist_name;
        const std::string xylabel_par = "hpar_xy_" + hist_name;

        const std::string xyzlabel = "hnu_xyz_" + hist_name;
        const std::string xyzlabel_par = "hpar_xyz_" + hist_name;

        hnu_E = std::make_unique<TH1D>(elabel.c_str(), elabel.c_str(), nbinsE, elow, eup);
        hnu_E_pipm = std::make_unique<TH1D>(elabel_pipm.c_str(), elabel_pipm.c_str(), nbinsE, elow, eup);
        hnu_E_kpm = std::make_unique<TH1D>(elabel_kpm.c_str(), elabel_kpm.c_str(), nbinsE, elow, eup);
        hnu_E_k0l = std::make_unique<TH1D>(elabel_k0l.c_str(), elabel_k0l.c_str(), nbinsE, elow, eup);

        hnu_theta_pipm = std::make_unique<TH1D>(thetalabel_pipm.c_str(), thetalabel_pipm.c_str(), nbinsTh, thlow, thup);
        hnu_theta_kpm = std::make_unique<TH1D>(thetalabel_kpm.c_str(), thetalabel_kpm.c_str(), nbinsTh, thlow, thup);
        hnu_theta_k0l = std::make_unique<TH1D>(thetalabel_k0l.c_str(), thetalabel_k0l.c_str(), nbinsTh, thlow, thup);

        hnu_thetaE_pipm =
            std::make_unique<TH2D>(thetaElabel_pipm.c_str(), thetaElabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
        hnu_thetaE_kpm =
            std::make_unique<TH2D>(thetaElabel_kpm.c_str(), thetaElabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
        hnu_thetaE_k0l =
            std::make_unique<TH2D>(thetaElabel_k0l.c_str(), thetaElabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);

        hnu_zE_pipm = std::make_unique<TH2D>(zElabel_pipm.c_str(), zElabel_pipm.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
        hnu_zE_kpm = std::make_unique<TH2D>(zElabel_kpm.c_str(), zElabel_kpm.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
        hnu_zE_k0l = std::make_unique<TH2D>(zElabel_k0l.c_str(), zElabel_k0l.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);

        hnu_yE_pipm = std::make_unique<TH2D>(yElabel_pipm.c_str(), yElabel_pipm.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
        hnu_yE_kpm = std::make_unique<TH2D>(yElabel_kpm.c_str(), yElabel_kpm.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
        hnu_yE_k0l = std::make_unique<TH2D>(yElabel_k0l.c_str(), yElabel_k0l.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);

        hnu_zy = std::make_unique<TH2D>(zylabel.c_str(), zylabel.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hnu_zx = std::make_unique<TH2D>(zxlabel.c_str(), zxlabel.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hnu_xy = std::make_unique<TH2D>(xylabel.c_str(), xylabel.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

        hnu_zy_pipm = std::make_unique<TH2D>(zylabel_pipm.c_str(), zylabel_pipm.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hnu_zx_pipm = std::make_unique<TH2D>(zxlabel_pipm.c_str(), zxlabel_pipm.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hnu_xy_pipm = std::make_unique<TH2D>(xylabel_pipm.c_str(), xylabel_pipm.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

        hnu_zy_kpm = std::make_unique<TH2D>(zylabel_kpm.c_str(), zylabel_kpm.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hnu_zx_kpm = std::make_unique<TH2D>(zxlabel_kpm.c_str(), zxlabel_kpm.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hnu_xy_kpm = std::make_unique<TH2D>(xylabel_kpm.c_str(), xylabel_kpm.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

        hnu_zy_k0l = std::make_unique<TH2D>(zylabel_k0l.c_str(), zylabel_k0l.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hnu_zx_k0l = std::make_unique<TH2D>(zxlabel_k0l.c_str(), zxlabel_k0l.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hnu_xy_k0l = std::make_unique<TH2D>(xylabel_k0l.c_str(), xylabel_k0l.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

        hpar_zy = std::make_unique<TH2D>(zylabel_par.c_str(), zylabel_par.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hpar_zx = std::make_unique<TH2D>(zxlabel_par.c_str(), zxlabel_par.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hpar_xy = std::make_unique<TH2D>(xylabel_par.c_str(), xylabel_par.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

        hnu_xyz = std::make_unique<TH3D>(xyzlabel.c_str(), xyzlabel.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
        hpar_xyz = std::make_unique<TH3D>(xyzlabel_par.c_str(), xyzlabel_par.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
    };

    void WriteHistograms();

    std::unique_ptr<TH1D> hnu_E;
    std::unique_ptr<TH1D> hnu_E_pipm;
    std::unique_ptr<TH1D> hnu_E_kpm;
    std::unique_ptr<TH1D> hnu_E_k0l;

    std::unique_ptr<TH1D> hnu_theta_pipm;
    std::unique_ptr<TH1D> hnu_theta_kpm;
    std::unique_ptr<TH1D> hnu_theta_k0l;

    std::unique_ptr<TH2D> hnu_thetaE_pipm;
    std::unique_ptr<TH2D> hnu_thetaE_kpm;
    std::unique_ptr<TH2D> hnu_thetaE_k0l;

    std::unique_ptr<TH2D> hnu_zE_pipm;
    std::unique_ptr<TH2D> hnu_zE_kpm;
    std::unique_ptr<TH2D> hnu_zE_k0l;

    std::unique_ptr<TH2D> hnu_yE_pipm;
    std::unique_ptr<TH2D> hnu_yE_kpm;
    std::unique_ptr<TH2D> hnu_yE_k0l;

    std::unique_ptr<TH2D> hnu_xy;
    std::unique_ptr<TH2D> hnu_zx;
    std::unique_ptr<TH2D> hnu_zy;

    std::unique_ptr<TH2D> hnu_xy_pipm;
    std::unique_ptr<TH2D> hnu_zx_pipm;
    std::unique_ptr<TH2D> hnu_zy_pipm;

    std::unique_ptr<TH2D> hnu_xy_kpm;
    std::unique_ptr<TH2D> hnu_zx_kpm;
    std::unique_ptr<TH2D> hnu_zy_kpm;

    std::unique_ptr<TH2D> hnu_xy_k0l;
    std::unique_ptr<TH2D> hnu_zx_k0l;
    std::unique_ptr<TH2D> hnu_zy_k0l;

    std::unique_ptr<TH2D> hpar_xy;
    std::unique_ptr<TH2D> hpar_zx;
    std::unique_ptr<TH2D> hpar_zy;

    std::unique_ptr<TH3D> hnu_xyz;
    std::unique_ptr<TH3D> hpar_xyz;
};

void Spectra::WriteHistograms()
{
    hnu_E->Write();
    hnu_E_pipm->Write();
    hnu_E_kpm->Write();
    hnu_E_k0l->Write();

    hnu_theta_pipm->Write();
    hnu_theta_kpm->Write();
    hnu_theta_k0l->Write();

    hnu_thetaE_pipm->Write();
    hnu_thetaE_kpm->Write();
    hnu_thetaE_k0l->Write();

    hnu_thetaE_pipm->Write();
    hnu_thetaE_kpm->Write();
    hnu_thetaE_k0l->Write();

    hnu_zE_pipm->Write();
    hnu_zE_kpm->Write();
    hnu_zE_k0l->Write();

    hnu_yE_pipm->Write();
    hnu_yE_kpm->Write();
    hnu_yE_k0l->Write();

    hnu_xy->Write();
    hnu_zx->Write();
    hnu_zy->Write();

    hnu_xy_pipm->Write();
    hnu_zx_pipm->Write();
    hnu_zy_pipm->Write();

    hnu_xy_kpm->Write();
    hnu_zx_kpm->Write();
    hnu_zy_kpm->Write();

    hnu_xy_k0l->Write();
    hnu_zx_k0l->Write();
    hnu_zy_k0l->Write();

    hpar_xy->Write();
    hpar_zx->Write();
    hpar_zy->Write();

    hnu_xyz->Write();
    hpar_xyz->Write();
}

void runEventLoop(TChain &chain, const bsim::Dk2Nu* dk2nu, Spectra &spec)
{
    chain.SetBranchAddress("dk2nu", &dk2nu);
    const std::array<double, 3> ICARUSCoords{450.37, 7991.98, 79512.66};
    const unsigned int n_entries = chain.GetEntries();

    double pct = 0.;

    static const double RAD2DEG = 180. / 3.1416;

    for (std::size_t entry = 0; entry < n_entries; ++entry)
    {
        pct = 100. * entry / n_entries;

        chain.GetEntry(entry);

        // only checking numu + numubar right now
        if (std::abs(dk2nu->decay.ntype) != 14)
        {
            continue;
        }

        auto const ptype = dk2nu->decay.ptype;
        const bool is_pipm = std::abs(ptype) == 211;
        const bool is_kpm = std::abs(ptype) == 321;
        const bool is_k0l = ptype == 130;

        double nu_energy = -9999.; // GeV
        double theta_par = -9999.; // rad

        const double wght = GetWeight(*dk2nu, ICARUSCoords, nu_energy, theta_par);

        theta_par *= RAD2DEG;

        if ((entry == 0) || (entry % 10000 == 0))
            std::cout << "FILLING ENTRY " << entry << " / " << n_entries << " (" << pct << "%)"
                      << " WITH WEIGHT = " << wght << " NU ENERGY = " << nu_energy << " GeV\n";

        const double nu_vx = dk2nu->decay.vx;
        const double nu_vy = dk2nu->decay.vy;
        const double nu_vz = dk2nu->decay.vz;

        const double par_vx = dk2nu->ppvx;
        const double par_vy = dk2nu->ppvy;
        const double par_vz = dk2nu->ppvz;

        spec.hnu_E->Fill(nu_energy, wght);

        spec.hnu_zy->Fill(nu_vz, nu_vy, wght);
        spec.hnu_zx->Fill(nu_vz, nu_vx, wght);
        spec.hnu_xy->Fill(nu_vx, nu_vy, wght);
        spec.hnu_xyz->Fill(nu_vx, nu_vy, nu_vz, wght);

        spec.hpar_zy->Fill(par_vz, par_vy, wght);
        spec.hpar_zx->Fill(par_vz, par_vx, wght);
        spec.hpar_xy->Fill(par_vx, par_vy, wght);
        spec.hpar_xyz->Fill(par_vx, par_vy, par_vz, wght);

        if (is_pipm)
        {
            spec.hnu_E_pipm->Fill(nu_energy, wght);
            spec.hnu_zy_pipm->Fill(nu_vz, nu_vy, wght);
            spec.hnu_zx_pipm->Fill(nu_vz, nu_vx, wght);
            spec.hnu_xy_pipm->Fill(nu_vx, nu_vy, wght);
            spec.hnu_theta_pipm->Fill(theta_par, wght);
            spec.hnu_thetaE_pipm->Fill(theta_par, nu_energy, wght);
            spec.hnu_yE_pipm->Fill(nu_vy, nu_energy, wght);
            spec.hnu_zE_pipm->Fill(nu_vz, nu_energy, wght);
        }
        else if (is_kpm)
        {
            spec.hnu_E_kpm->Fill(nu_energy, wght);
            spec.hnu_zy_kpm->Fill(nu_vz, nu_vy, wght);
            spec.hnu_zx_kpm->Fill(nu_vz, nu_vx, wght);
            spec.hnu_xy_kpm->Fill(nu_vx, nu_vy, wght);
            spec.hnu_theta_kpm->Fill(theta_par, wght);
            spec.hnu_thetaE_kpm->Fill(theta_par, nu_energy, wght);
            spec.hnu_yE_kpm->Fill(nu_vy, nu_energy, wght);
            spec.hnu_zE_kpm->Fill(nu_vz, nu_energy, wght);
        }
        else if (is_k0l)
        {
            spec.hnu_E_k0l->Fill(nu_energy, wght);
            spec.hnu_zy_k0l->Fill(nu_vz, nu_vy, wght);
            spec.hnu_zx_k0l->Fill(nu_vz, nu_vx, wght);
            spec.hnu_xy_k0l->Fill(nu_vx, nu_vy, wght);
            spec.hnu_theta_k0l->Fill(theta_par, wght);
            spec.hnu_thetaE_k0l->Fill(theta_par, nu_energy, wght);
            spec.hnu_yE_k0l->Fill(nu_vy, nu_energy, wght);
            spec.hnu_zE_k0l->Fill(nu_vz, nu_energy, wght);
        }
    } // Event Loop
} // fillSpectra

int main()
{
    const double pot_per_file = 500000.;

    auto const files_without_blocks = getFileList("files_without_blocks.txt");
    auto const files_with_blocks = getFileList("kaon_xsec_files.txt");
    // auto const files_with_blocks = getFileList("files_with_blocks.txt");

    auto pChain_no_blocks = std::make_unique<TChain>("dk2nuTree");
    auto pChain_blocks = std::make_unique<TChain>("dk2nuTree");

    for (auto const &f : files_without_blocks)
    {
        pChain_no_blocks->Add(f.c_str());
    }

    for (auto const &f : files_with_blocks)
    {
        pChain_blocks->Add(f.c_str());
    }

    auto pDk2nu_no_blocks = std::make_unique<bsim::Dk2Nu>();
    auto pDk2nu_blocks = std::make_unique<bsim::Dk2Nu>();

    Spectra spec_no_blocks("no_blocks");
    Spectra spec_blocks("blocks");

    // g3Chase = false
    runEventLoop(*pChain_no_blocks, pDk2nu_no_blocks.get(), spec_no_blocks);

    // g3Chase = true
    runEventLoop(*pChain_blocks, pDk2nu_blocks.get(), spec_blocks);

    const double total_pot_blocks = pot_per_file * (double)files_with_blocks.size();
    const double total_pot_no_blocks = pot_per_file * (double)files_without_blocks.size();

    auto hpot_no_blocks = std::make_unique<TH1D>("hpot_no_blocks", "hpot_no_blocks", 1, 0, 1);
    auto hpot_blocks = std::make_unique<TH1D>("hpot_blocks", "hpot_blocks", 1, 0, 1);

    hpot_no_blocks->SetBinContent(1, total_pot_no_blocks);
    hpot_blocks->SetBinContent(1, total_pot_blocks);

    std::cout << "\nSaving hists to out.root...\n";

    TFile fOut("out.root", "RECREATE");

    spec_no_blocks.WriteHistograms();
    spec_blocks.WriteHistograms();

    hpot_no_blocks->Write();
    hpot_blocks->Write();

    fOut.Close();

    std::cout << "Done." << std::endl;

    return 0;
}
