#ifndef PDGID_H
#define PDGID_H

#include <unordered_map>

namespace pdg {
  enum PDGID {
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
    NULLTARGET = 1000000000, // Decay process ID in geant4.10+
    ZERO = 0,                // Decay process ID in old versions
  };

  static const std::unordered_map<PDGID, double> pdgid2Mass = {{PIP, 0.13957039},
                                                               {PIM, .13957039},
                                                               {KP, 0.493677},
                                                               {KM, 0.493677},
                                                               {K0, 0.497611},
                                                               {K0L, 0.497611},
                                                               {K0S, 0.497611},
                                                               {MUM, 0.1056583755},
                                                               {MUP, 0.1056583755},
                                                               {TAUM, 1.77686},
                                                               {TAUP, 1.77686}};
} //namespace pdg

#endif // PDGID_H
