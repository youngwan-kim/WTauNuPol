#include <TauPolarization/TauPolarimeter.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace TauPol {

int GenDecayModeFromCounts(int nChargedHadrons, int nPi0) {
    if (nChargedHadrons == 1 && nPi0 >= 0 && nPi0 <= 2)
        return nPi0;
    if (nChargedHadrons == 3 && (nPi0 == 0 || nPi0 == 1))
        return 10 + nPi0;
    return kDMOther;
}

double PdgMass(int pdgId, double storedMass) {
    if (storedMass > 1e-3)
        return storedMass;
    switch (std::abs(pdgId)) {
    case 11: return 0.000510999;
    case 13: return 0.1056584;
    case 15: return kTauMass;
    case 111: return kPiNeutralMass;
    case 211: return kPiChargedMass;
    case 321: return 0.493677;
    case 130:
    case 310:
    case 311: return 0.497611;
    case 221: return 0.547862;
    default: return storedMass;
    }
}

namespace {

void FillLabFractions(const TLorentzVector &tau, const TLorentzVector &charged, const TLorentzVector &visible, PolarimeterResult &out) {
    if (tau.E() > 0.)
        out.x = static_cast<float>(visible.E() / tau.E());
    if (visible.E() > 0. && charged.E() > 0.)
        out.upsilon = static_cast<float>(2. * charged.E() / visible.E() - 1.);
}

} // namespace

PolarimeterResult ComputePolarimeter(const TLorentzVector &tau, const TLorentzVector &charged, const TLorentzVector &neutral, int genDM) {
    PolarimeterResult out;
    out.genDM = genDM;
    const TLorentzVector visible = charged + neutral;
    FillLabFractions(tau, charged, visible, out);
    if (genDM != kDM0 && genDM != kDM1)
        return out;
    if (tau.P() <= 0. || charged.E() <= 0.)
        return out;

    // helicity frame
    const TVector3 boost = -tau.BoostVector();
    const TVector3 zhat = tau.Vect().Unit();
    TLorentzVector P = tau, C = charged, Nt = neutral;
    P.Boost(boost);
    C.Boost(boost);
    Nt.Boost(boost);
    const TLorentzVector Q = C + Nt;
    const TLorentzVector N = P - Q; 

    TLorentzVector J;
    if (genDM == kDM0) {
        J = C;
    } else {
        if (Nt.E() <= 0.)
            return out;
        const TLorentzVector q = C - Nt;
        const double Q2 = Q.M2();
        if (Q2 <= 0.)
            return out;
        J = q - (q.Dot(Q) / Q2) * Q;
    }

    const double JN = J.Dot(N);
    const double JP = J.Dot(P);
    const double JJ = J.Dot(J);
    const double NP = N.Dot(P);
    const double U = 2. * JN * JP - JJ * NP;
    const double m = P.M();
    const TVector3 S = m * (2. * JN * J.Vect() - JJ * N.Vect());
    if (!(U > 0.))
        return out;

    double omega = S.Dot(zhat) / U;
    
    out.clamped = std::abs(omega) > 1. + 1e-6;
    omega = std::clamp(omega, -1., 1.);

    if (genDM == kDM0) {
        // omega = cos(theta*) exactly. 
        // evaluated it from the lab energy fraction,
        //   x m_tau = E* + beta p* cos(theta*),
        // which needs no angles. NanoAOD stores GenPart pt/eta/phi with a
        // reduced mantissa, so the rest-frame angle built from directions is
        // smeared by O(gamma x 1e-3) 
        const double mc2 = C.M2() > 0. ? C.M2() : 0.;
        const double eStar = (m * m + mc2) / (2. * m);
        const double pStar = (m * m - mc2) / (2. * m);
        const double beta = tau.Beta();
        if (!(pStar > 0.) || !(beta > 0.))
            return out;
        omega = (charged.E() / tau.E() * m - eStar) / (beta * pStar);
        out.clamped = std::abs(omega) > 1. + 1e-6;
        omega = std::clamp(omega, -1., 1.);
    }

    out.valid = true;
    out.u = U;
    out.omega = static_cast<float>(omega);
    if (genDM == kDM0)
        out.cosThetaStar = out.omega;
    else if (Q.P() > 0.)
        out.cosThetaStar = static_cast<float>(Q.Vect().Unit().Dot(zhat));
    if (genDM == kDM1 && Q.P() > 0.) {
        TLorentzVector cInQ = C;
        cInQ.Boost(-Q.BoostVector());
        if (cInQ.P() > 0.)
            out.cosPsi = static_cast<float>(cInQ.Vect().Unit().Dot(Q.Vect().Unit()));
    }
    return out;
}

double RawTemplateWeight(float omega, float pTarget, float pGen) {
    if (pTarget == pGen)
        return 1.; // also at 1 + pGen omega = 0 (omega = +-1 exactly), where the ratio is 0/0
    const double den = 1. + static_cast<double>(pGen) * omega;
    if (!(den > 0.))
        return std::numeric_limits<double>::infinity();
    return (1. + static_cast<double>(pTarget) * omega) / den;
}

float TemplateWeight(float omega, float pTarget, float pGen, float cap) {
    return static_cast<float>(std::min<double>(RawTemplateWeight(omega, pTarget, pGen), cap));
}

float UpsilonFromRatio(float leadTkPtOverTauPt) { return 2.f * leadTkPtOverTauPt - 1.f; }

float ShiftedRatio(float ratio, float miniAODPt, float shiftedPt) {
    if (!(shiftedPt > 0.f))
        return ratio;
    return ratio * miniAODPt / shiftedPt;
}

double ShapeFisherInformation(const std::vector<double> &nPlus, const std::vector<double> &nMinus, double p) {
    const std::size_t nb = std::min(nPlus.size(), nMinus.size());
    double A = 0., A1 = 0.;
    for (std::size_t b = 0; b < nb; ++b) {
        const double f0 = 0.5 * (nPlus[b] + nMinus[b]);
        const double f1 = 0.5 * (nPlus[b] - nMinus[b]);
        if (f0 + p * f1 <= 0.)
            continue;
        A += f0 + p * f1;
        A1 += f1;
    }
    if (!(A > 0.))
        return 0.;
    double info = 0.;
    for (std::size_t b = 0; b < nb; ++b) {
        const double f0 = 0.5 * (nPlus[b] + nMinus[b]);
        const double f1 = 0.5 * (nPlus[b] - nMinus[b]);
        const double fp = f0 + p * f1;
        if (fp <= 0.)
            continue;
        const double d = f1 / fp - A1 / A;
        info += (fp / A) * d * d;
    }
    return info;
}

} // namespace TauPol
