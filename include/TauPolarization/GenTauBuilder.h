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

enum class TauOrigin : int { W = 0, TopW = 1 };

struct GenTau {
    int index = -1;           // last-copy tau in GenPart
    int charge = 0;           // +-1
    TauOrigin origin = TauOrigin::W;
    int ancestorPdgId = 0;    // first non tau ancestor (+-24)
    TLorentzVector tau, charged, neutral, visible;
    int chargedPdgId = 0;     // leading charged hadron 
    int nChargedHadrons = 0, nPi0 = 0, nPhotons = 0, nOtherNeutral = 0;
    bool leptonic = false;
    int genDM = kDMOther;     // from the stored daughters
    bool hasNu = false;
    TLorentzVector nu;
};

constexpr unsigned short kIsLastCopy = 1u << 13;

// all last copy taus whose first non tau ancestor is a W, taus from a W whose 
// first non-W ancestor is a top are tagged as TopW.
// taus from other sources are counted as nOtherOrigin 
std::vector<GenTau> BuildGenTausFromW(const GenRecord &gen, int *nOtherOrigin = nullptr);

} // namespace TauPol

#endif
