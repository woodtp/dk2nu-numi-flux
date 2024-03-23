#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStopwatch.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
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

enum PDGID
{
    PIP = 211,
    PIM = -211,
    KP = 321,
    KM = -321,
    K0 = 311,
    K0L = 130,
    K0S = 310,
    MUM = 13,
    MUP = -13,
    TAUM = 15,
    TAUP = -15,
    NUE = 12,
    NUEBAR = -12,
    NUMU = 14,
    NUMUBAR = -14,
    NUTAU = 16,
    NUTAUBAR = -16,
    PROTON = 2212,
    NEUTRON = 2112,
    CARBON = 1000060120,
    ALUMINUM = 1000130270,
    IRON = 1000260560,
    NULLTARGET = 1000000000,
    ZERO = 0,
};

static const std::map<PDGID, double> pdgid2Mass = {{PIP, 0.13957039}, {PIM, .13957039}, {KP, 0.493677}, {KM, 0.493677}, {K0, 0.497611},
    {K0L, 0.497611}, {K0S, 0.497611}, {MUM, 0.1056583755}, {MUP, 0.1056583755}, {TAUM, 1.77686}, {TAUP, 1.77686}};

double GetWeight(const bsim::Dk2Nu &dk2nu, const std::array<double, 3> &detCoords, double &nu_energy)
{
    static constexpr double rdet = 100.0; //in cm

    const auto parent_ptype = static_cast<PDGID>(dk2nu.decay.ptype);

    if (pdgid2Mass.find(parent_ptype) == pdgid2Mass.end())
    {
        std::cout << "GetWeight - Wrong parent type!! " << parent_ptype << '\n';
        return -9999.;
    }

    const double parent_mass = pdgid2Mass.at(parent_ptype);

    const double pdPx = dk2nu.decay.pdpx;
    const double pdPy = dk2nu.decay.pdpy;
    const double pdPz = dk2nu.decay.pdpz;

    const double parent_energy = std::sqrt((pdPx * pdPx) + (pdPy * pdPy) + (pdPz * pdPz) + (parent_mass * parent_mass));
    const double gamma = parent_energy / parent_mass;
    const double gamma_sqr = gamma * gamma;
    const double beta = std::sqrt((gamma_sqr - 1.) / gamma_sqr);

    const double parentP = std::sqrt((pdPx * pdPx) + (pdPy * pdPy) + (pdPz * pdPz));

    const double rrx = detCoords[0] - dk2nu.decay.vx;
    const double rry = detCoords[1] - dk2nu.decay.vy;
    const double rrz = detCoords[2] - dk2nu.decay.vz;
    const double rr = std::sqrt(rrx * rrx + rry * rry + rrz * rrz);

    double costh_parent = ((pdPx * rrx) + (pdPy * rry) + (pdPz * rrz)) / (parentP * rr);

    if (costh_parent > 1.)
    {
        costh_parent = 1.;
    }
    else if (costh_parent < -1.)
    {
        costh_parent = -1.;
    }
    const double emrat = 1. / (gamma * (1. - beta * costh_parent));
    const double angle_detector = (rdet * rdet) / (rr * rr) / 4.;

    double wght = (angle_detector * emrat * emrat * dk2nu.decay.nimpwt) / 3.1416;

    const double enuzr = dk2nu.decay.necm;

    nu_energy = emrat * enuzr;

    //done for all except polarized muon
    // in which case need to modify weight
    if (std::abs(parent_ptype) == PDGID::MUM)
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
            if (nu_type == PDGID::NUE || nu_type == PDGID::NUEBAR)
            {
                wt_ratio = 1. - costh;
            }
            else if (nu_type == PDGID::NUMU || nu_type == PDGID::NUMUBAR)
            {
                const double mumass = pdgid2Mass.at(MUM);
                double xnu = 2. * enuzr / mumass;
                wt_ratio = ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
            }
            else if (nu_type == PDGID::NUTAU || nu_type == PDGID::NUTAUBAR)
            {
                const double taumass = pdgid2Mass.at(PDGID::TAUM);
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
    Spectra(const std::string &id, const int nu_pdg = -1)
    {
        pdg2Index = {{PROTON, 0.5}, {NEUTRON, 1.5}, {PIP, 2.5}, {PIM, 3.5}, {KP, 4.5}, {KM, 5.5}, {K0L, 6.5}, {CARBON, 0.5}, {ALUMINUM, 1.5}, {IRON, 2.5}, {NULLTARGET, 3.5}, {ZERO, 4.5}};

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
        if (nu_pdg < 0 && nu_pdg != -1)
        {
            nu_suff += "bar";
        }

        const std::string suffix = id + "_" + nu_suff;

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
        hancestorInteractions = TH2F(ints_label.c_str(), ints_label.c_str(), 6, 0., 6., 8, 0., 8.);

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

        hnu_thetaE = TH2D(thetaElabel.c_str(), thetaElabel.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
        hnu_thetaE_pipm = TH2D(thetaElabel_pipm.c_str(), thetaElabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
        hnu_thetaE_kpm = TH2D(thetaElabel_kpm.c_str(), thetaElabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
        hnu_thetaE_k0l = TH2D(thetaElabel_k0l.c_str(), thetaElabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);
        hnu_thetaE_mu = TH2D(thetaElabel_mu.c_str(), thetaElabel_mu.c_str(), nbinsTh, thlow, thup, nbinsE, elow, eup);

        hnu_thetaY = TH2D(thetaYlabel.c_str(), thetaYlabel.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
        hnu_thetaY_pipm = TH2D(thetaYlabel_pipm.c_str(), thetaYlabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
        hnu_thetaY_kpm = TH2D(thetaYlabel_kpm.c_str(), thetaYlabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
        hnu_thetaY_k0l = TH2D(thetaYlabel_k0l.c_str(), thetaYlabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);
        hnu_thetaY_mu = TH2D(thetaYlabel_mu.c_str(), thetaYlabel_mu.c_str(), nbinsTh, thlow, thup, nbinsy, ylow, yup);

        hnu_thetaZ = TH2D(thetaZlabel.c_str(), thetaZlabel.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
        hnu_thetaZ_pipm = TH2D(thetaZlabel_pipm.c_str(), thetaZlabel_pipm.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
        hnu_thetaZ_kpm = TH2D(thetaZlabel_kpm.c_str(), thetaZlabel_kpm.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
        hnu_thetaZ_k0l = TH2D(thetaZlabel_k0l.c_str(), thetaZlabel_k0l.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);
        hnu_thetaZ_mu = TH2D(thetaZlabel_mu.c_str(), thetaZlabel_mu.c_str(), nbinsTh, thlow, thup, nbinsz, zlow, zup);

        hnu_zE = TH2D(zElabel.c_str(), zElabel.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
        hnu_zE_pipm = TH2D(zElabel_pipm.c_str(), zElabel_pipm.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
        hnu_zE_kpm = TH2D(zElabel_kpm.c_str(), zElabel_kpm.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
        hnu_zE_k0l = TH2D(zElabel_k0l.c_str(), zElabel_k0l.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);
        hnu_zE_mu = TH2D(zElabel_mu.c_str(), zElabel_mu.c_str(), nbinsz, zlow, zup, nbinsE, elow, eup);

        hnu_yE = TH2D(yElabel.c_str(), yElabel.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
        hnu_yE_pipm = TH2D(yElabel_pipm.c_str(), yElabel_pipm.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
        hnu_yE_kpm = TH2D(yElabel_kpm.c_str(), yElabel_kpm.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
        hnu_yE_k0l = TH2D(yElabel_k0l.c_str(), yElabel_k0l.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);
        hnu_yE_mu = TH2D(yElabel_mu.c_str(), yElabel_mu.c_str(), nbinsy, ylow, yup, nbinsE, elow, eup);

        hnu_zy = TH2D(zylabel.c_str(), zylabel.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hnu_zx = TH2D(zxlabel.c_str(), zxlabel.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hnu_xy = TH2D(xylabel.c_str(), xylabel.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

        hnu_zy_pipm = TH2D(zylabel_pipm.c_str(), zylabel_pipm.c_str(), nbinsz, zlow, zup, nbinsy, ylow, yup);
        hnu_zx_pipm = TH2D(zxlabel_pipm.c_str(), zxlabel_pipm.c_str(), nbinsz, zlow, zup, nbinsx, xlow, xup);
        hnu_xy_pipm = TH2D(xylabel_pipm.c_str(), xylabel_pipm.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);

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

        hnu_xyz = TH3D(xyzlabel.c_str(), xyzlabel.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
        hpar_xyz = TH3D(xyzlabel_par.c_str(), xyzlabel_par.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
    };

    void FillSpectra(const bsim::Dk2Nu &dk2nu, const double wght, const double nu_energy, const double theta_par);
    void WriteHistograms();

private:
    long unsigned int nNeutrinos = 0;
    long unsigned int nInteractions = 0;


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

    std::map<PDGID, float> pdg2Index;
    std::vector<std::pair<PDGID, PDGID>> IdentifyPrecursors(const bsim::Dk2Nu &dk2nu);
};

std::vector<std::pair<PDGID, PDGID>> Spectra::IdentifyPrecursors(const bsim::Dk2Nu &dk2nu)
{
    std::vector<std::pair<PDGID, PDGID>> precursors;
    precursors.reserve(dk2nu.ancestor.size());
    // for(auto const& ancestor : dk2nu.ancestor)
    for(auto i = dk2nu.ancestor.begin(); i != dk2nu.ancestor.end() - 1; ++i)
    {
        const auto parent_ptype = static_cast<PDGID>(i->pdg);
        const auto target_type = static_cast<PDGID>(i->nucleus);
        precursors.emplace_back(std::make_pair(parent_ptype, target_type));
    }
    return precursors;
}

void Spectra::FillSpectra(const bsim::Dk2Nu &dk2nu, const double wght, const double nu_energy, const double theta_par)
{

    const double nu_vx = dk2nu.decay.vx;
    const double nu_vy = dk2nu.decay.vy;
    const double nu_vz = dk2nu.decay.vz;

    const double par_vx = dk2nu.ppvx;
    const double par_vy = dk2nu.ppvy;
    const double par_vz = dk2nu.ppvz;

    if ( nu_energy > 0.4 )
    {
        nNeutrinos++;

        auto const& precursors = IdentifyPrecursors(dk2nu);

        for (auto const& precursor : precursors)
        {
            auto const &target = precursor.second;
            auto const &proj = precursor.first;

            // NULLTARGET occurs when there is a decay in flight or at rest, or when the neutrino is produced
            // if (target == NULLTARGET || (proj != PROTON && target == 0))
            // {
            //     continue;
            // }

            nInteractions++;

            auto const projIdx = pdg2Index.find(proj) != pdg2Index.end() ? pdg2Index.at(proj) : 7.5;
            auto const targetIdx = pdg2Index.find(target) != pdg2Index.end() ? pdg2Index.at(target) : 5.5;
            hancestorInteractions.Fill(targetIdx, projIdx, wght);
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

    auto const ptype = std::abs(dk2nu.decay.ptype);

    if (ptype == PDGID::PIP)
    {
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
    else if (ptype == PDGID::KP)
    {
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
    else if (ptype == PDGID::K0L)
    {
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
    else if (ptype == PDGID::MUM)
    {
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
    hnNeutrinos.SetBinContent(1, nNeutrinos);
    hnInteractions.SetBinContent(1, nInteractions);

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

double mag(const std::array<double, 3> &u)
{
    return std::sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]));
}

double dot(const std::array<double, 3> &u, const std::array<double, 3> &v)
{
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

std::array<double, 3> cross(const std::array<double, 3> &u, const std::array<double, 3> &v)
{
    const std::array<double, 3> res = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    return res;
}

std::array<double, 3> project(const std::array<double, 3> &u, const std::array<double, 3> &v)
{
    const double vMag = mag(v);
    const double costh = dot(u, v) / (vMag * vMag);
    const std::array<double, 3> res = {u[0] - costh * v[0], u[1] - costh * v[1], u[2] - costh * v[2]};
    return res;
}

double calculateTheta(const std::array<double, 3> &u, const std::array<double, 3> &v)
{
    const double uMag = mag(u);
    const double vMag = mag(v);
    const double costh = dot(u, v) / (uMag * vMag);
    if (costh > 1.)
    {
        return std::acos(1.);
    }
    else if (costh < -1.)
    {
        return std::acos(-1.);
    }
    return std::acos(costh); // [rad]
}

void runEventLoop(TChain &chain, const bsim::Dk2Nu *dk2nu, std::map<int, Spectra *> &spectra)
{
    chain.SetBranchAddress("dk2nu", &dk2nu);

    static const std::array<double, 3> ICARUSCoords = {450.37, 7991.98, 79512.66};
    static const std::array<double, 3> NuMIAxis = {0., 0., 1.};

    static const auto n = cross(ICARUSCoords, NuMIAxis);
    static constexpr double magN = 1.;

    const unsigned int n_entries = chain.GetEntries();

    double pct = 0.;

    static constexpr double RAD2DEG = 180. / 3.1416;

    for (std::size_t entry = 0; entry < n_entries; ++entry)
    {
        pct = 100. * entry / n_entries;

        chain.GetEntry(entry);
        auto &spec = spectra[dk2nu->decay.ntype];
        double nu_energy = -9999.; // GeV
        const std::array<double, 3> parentDecayP = {dk2nu->decay.pdpx, dk2nu->decay.pdpy, dk2nu->decay.pdpz};
        const auto parentDecayPProj = project(parentDecayP, n);
        const auto magPProj = mag(parentDecayPProj);

        const double costheta_par = dot(parentDecayPProj, NuMIAxis) / (magPProj * magN);

        double theta_par = RAD2DEG * std::acos(costheta_par);

        // check if parentDecayPProj is between NuMIAxis and ICARUSCoords on the plane formed by NuMIAxis and ICARUSCoords
        // and negate theta if it is not.
        const auto crossProd = cross(parentDecayPProj, NuMIAxis);
        const double dotProd = dot(crossProd, n);

        if (dotProd < 0.)
        {
            theta_par *= -1.;
        }

        const double wght = GetWeight(*dk2nu, ICARUSCoords, nu_energy);

        spec->FillSpectra(*dk2nu, wght, nu_energy, theta_par);

        if ((entry == 0) || (entry % 10000 == 0))
        {
            std::cout << "FILLING ENTRY " << entry << " / " << n_entries << " (" << pct << "%)"
                      << " WITH WEIGHT = " << wght << '\n';
        }

    } // Event Loop
} // runEventLoop

int main()
{
    TStopwatch sw;
    sw.Start();

    ROOT::EnableImplicitMT(1000);

    static constexpr double pot_per_file = 500000.;

    // // auto const files_without_blocks = getFileList("files_without_blocks.txt");
    // // auto const files_with_blocks = getFileList("files_with_blocks.txt");
    // // auto const files_with_blocks_kaons = getFileList("kaon_xsec_files.txt");

    auto const files_without_blocks = getFileList("default_files.txt");
    auto const files_with_blocks = getFileList("g3Chase_files.txt");
    auto const files_with_blocks_kaons = getFileList("updated_g4_filelist_50M_POT.txt");

    auto pChain_no_blocks = std::make_unique<TChain>("dk2nuTree");
    auto pChain_blocks = std::make_unique<TChain>("dk2nuTree");
    auto pChain_blocks_kaons = std::make_unique<TChain>("dk2nuTree");

    for (auto const &f : files_without_blocks)
    {
        pChain_no_blocks->Add(f.c_str());
    }

    for (auto const &f : files_with_blocks)
    {
        pChain_blocks->Add(f.c_str());
    }

    for (auto const &f : files_with_blocks_kaons)
    {
        pChain_blocks_kaons->Add(f.c_str());
    }

    const bsim::Dk2Nu dk2nu_no_blocks;
    const bsim::Dk2Nu dk2nu_blocks;
    const bsim::Dk2Nu dk2nu_blocks_kaons;

    Spectra spec_fhc_numu_no_blocks("no_blocks", 14);
    Spectra spec_fhc_numubar_no_blocks("no_blocks", -14);
    Spectra spec_fhc_nue_no_blocks("no_blocks", 12);
    Spectra spec_fhc_nuebar_no_blocks("no_blocks", -12);

    Spectra spec_fhc_numu_blocks("blocks", 14);
    Spectra spec_fhc_numubar_blocks("blocks", -14);
    Spectra spec_fhc_nue_blocks("blocks", 12);
    Spectra spec_fhc_nuebar_blocks("blocks", -12);

    Spectra spec_fhc_numu_blocks_kaons("blocks_kaons", 14);
    Spectra spec_fhc_numubar_blocks_kaons("blocks_kaons", -14);
    Spectra spec_fhc_nue_blocks_kaons("blocks_kaons", 12);
    Spectra spec_fhc_nuebar_blocks_kaons("blocks_kaons", -12);

    std::map<int, Spectra *> spectra_no_blocks = {{14, &spec_fhc_numu_no_blocks}, {-14, &spec_fhc_numubar_no_blocks},
        {12, &spec_fhc_nue_no_blocks}, {-12, &spec_fhc_nuebar_no_blocks}};

    std::map<int, Spectra *> spectra_blocks = {
        {14, &spec_fhc_numu_blocks}, {-14, &spec_fhc_numubar_blocks}, {12, &spec_fhc_nue_blocks}, {-12, &spec_fhc_nuebar_blocks}};

    std::map<int, Spectra *> spectra_blocks_kaons = {{14, &spec_fhc_numu_blocks_kaons}, {-14, &spec_fhc_numubar_blocks_kaons},
        {12, &spec_fhc_nue_blocks_kaons}, {-12, &spec_fhc_nuebar_blocks_kaons}};

    const double total_pot_no_blocks = pot_per_file * (double)files_without_blocks.size();
    const double total_pot_blocks = pot_per_file * (double)files_with_blocks.size();
    const double total_pot_blocks_kaons = pot_per_file * (double)files_with_blocks_kaons.size();

    // g3Chase = false
    runEventLoop(*pChain_no_blocks, &dk2nu_no_blocks, spectra_no_blocks);

    // g3Chase = true
    runEventLoop(*pChain_blocks, &dk2nu_blocks, spectra_blocks);

    // g3Chase = true + kaon xsec bugfix
    runEventLoop(*pChain_blocks_kaons, &dk2nu_blocks_kaons, spectra_blocks_kaons);
    TH1D hpot_no_blocks("hpot_no_blocks", "hpot_no_blocks", 1, 0, 1);
    TH1D hpot_blocks("hpot_blocks", "hpot_blocks", 1, 0, 1);
    TH1D hpot_blocks_kaons("hpot_blocks_kaons", "hpot_blocks_kaons", 1, 0, 1);

    hpot_no_blocks.SetBinContent(1, total_pot_no_blocks);
    hpot_blocks.SetBinContent(1, total_pot_blocks);
    hpot_blocks_kaons.SetBinContent(1, total_pot_blocks_kaons);

    std::cout << "\nSaving hists to out.root...\n";
    TFile fOut("out.root", "RECREATE");

    for (auto &it : spectra_no_blocks)
    {
        it.second->WriteHistograms();
    }

    for (auto &it : spectra_blocks)
    {
        it.second->WriteHistograms();
    }

    for (auto &it : spectra_blocks_kaons)
    {
        it.second->WriteHistograms();
    }

    hpot_no_blocks.Write();
    hpot_blocks.Write();
    hpot_blocks_kaons.Write();

    fOut.Close();
    sw.Stop();
    std::cout << "Done in " << sw.RealTime() << " seconds." << std::endl;

    return 0;
}
