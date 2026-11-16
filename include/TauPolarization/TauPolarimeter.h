#ifndef TauPolarization_TauPolarimeter_h
#define TauPolarization_TauPolarimeter_h

#include <vector>
#include "TLorentzVector.h"

// per-tau polarimeter and template helpers

// conventions
// 1) metric (+,-,-,-), GeV. 
// 2) P_tau is quoted for tau-, the same formulas are applied to both charges 
//    (the charged pion of a tau+ decay is treated exactly as the pi- of a tau- decay), so omega measures P(tau-) = -P(tau+) for either charge.

namespace TauPol {

constexpr double kTauMass = 1.77686;
constexpr double kPiChargedMass = 0.13957039;
constexpr double kPiNeutralMass = 0.1349768;

enum GenDecayMode : int {
    kDMLeptonic = -1,
    kDM0 = 0,   // h nu
    kDM1 = 1,   // h pi0 nu
    kDM2 = 2,   // h 2pi0 nu
    kDM10 = 10, // 3h nu
    kDM11 = 11, // 3h pi0 nu
    kDMOther = 15,
};

// 5 * (nCharged - 1) + nPi0, with anything outside {0,1,2,10,11} is mapped to 15
int GenDecayModeFromCounts(int nChargedHadrons, int nPi0);

double PdgMass(int pdgId, double storedMass);

struct PolarimeterResult {
    bool valid = false;
    float omega = -9.f;         // S_z / U in [-1, 1]
    float cosThetaStar = -9.f;  // visible system (pi for DM0, rho for DM1) vs z, tau rest frame
    float cosPsi = -9.f;        // DM1: pi+- in the (pi+- + neutral) frame vs that system's direction
    float x = -9.f;             // E_vis / E_tau (lab)
    float upsilon = -9.f;       // 2 E_charged / E_vis - 1 (lab)
    int genDM = kDMOther;
    bool clamped = false;       // |omega| exceeded 1 (imprecise inputs) and was clamped
    double u = 0.;              // unpolarized density U (for toy reweighting); rate ~ U (1 + P omega)
};

// tau, charged, neutral are lab four-vectors
// neutral may be (0,0,0,0) for DM0. omega is computed for genDM 0 and 1 only, 
// x, upsilon and the angles are filled whenever the kinematics allow. 
// neutrino is defined as N = tau - (charged + neutral)
PolarimeterResult ComputePolarimeter(const TLorentzVector &tau, const TLorentzVector &charged, const TLorentzVector &neutral, int genDM);

// w(P) = (1 + P omega) / (1 + pGen omega), uncapped (+inf if the denominator is not positive).
double RawTemplateWeight(float omega, float pTarget, float pGen);
// min(RawTemplateWeight, cap). a weight is capped iff raw > cap.
float TemplateWeight(float omega, float pTarget, float pGen, float cap);

float UpsilonFromRatio(float leadTkPtOverTauPt);                    // 2R - 1
float ShiftedRatio(float ratio, float miniAODPt, float shiftedPt);  // R * Pt0 / Pt'

// shape-only Fisher information per event
// nPlus/nMinus are the weighted yields per bin of the P = +1 / -1 templates, acceptance included and not renormalized. 
// bins with F_b(P) <= 0 are skipped
double ShapeFisherInformation(const std::vector<double> &nPlus, const std::vector<double> &nMinus, double p);

} // namespace TauPol

#endif
