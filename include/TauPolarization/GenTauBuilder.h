#ifndef TauPolarization_GenTauBuilder_h
#define TauPolarization_GenTauBuilder_h

#include <cstddef>
#include <vector>
#include "TLorentzVector.h"
#include <TauPolarization/TauPolarimeter.h>

// reconstruct generator taus from W decays out of flat GenPart columns

namespace TauPol {

struct GenRecord {
    std::vector<float> pt, eta, phi, mass;
    std::vector<int> pdgId, mother;
    std::vector<unsigned short> statusFlags;
    std::size_t size() const { return pt.size(); }
    TLorentzVector P4(std::size_t i) const; // PDG mass for light particles
};

enum class TauOrigin : int { W = 0, TopW = 1, Z = 2 };

struct GenTau {
    int index = -1;           // last-copy tau in GenPart
    int charge = 0;           // +-1
    TauOrigin origin = TauOrigin::W;
    int ancestorPdgId = 0;    // first non tau ancestor (+-24, or 23 with acceptZ)
    TLorentzVector tau, charged, neutral, visible;
    int chargedPdgId = 0;     // leading charged hadron 
    int nChargedHadrons = 0, nPi0 = 0, nPhotons = 0, nOtherNeutral = 0;
    bool leptonic = false;
    int genDM = kDMOther;     // from the stored daughters
    bool hasNu = false;
    TLorentzVector nu;
};

constexpr unsigned short kIsLastCopy = 1u << 13;

std::vector<GenTau> BuildGenTausFromW(const GenRecord &gen, int *nOtherOrigin = nullptr);
std::vector<GenTau> BuildGenTaus(const GenRecord &gen, bool acceptZ, int *nOtherOrigin = nullptr);

} // namespace TauPol

#endif
