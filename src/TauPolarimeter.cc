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
    if (genDM == kDM0) {
        out.validE = true;
        out.omegaE = out.omega;
    } else {
        const double oe = OmegaEnergyDM1(tau.E(), m, visible.E(), charged.E(), visible.M(), charged.M() > 0. ? charged.M() : kPiChargedMass, neutral.M() > 0. ? neutral.M() : kPiNeutralMass);
        if (std::isfinite(oe)) {
            out.validE = true;
            out.omegaE = static_cast<float>(oe);
        }
    }
    if (genDM == kDM1 && Q.P() > 0.) {
        TLorentzVector cInQ = C;
        cInQ.Boost(-Q.BoostVector());
        if (cInQ.P() > 0.)
            out.cosPsi = static_cast<float>(cInQ.Vect().Unit().Dot(Q.Vect().Unit()));
    }
    return out;
}

namespace {
struct V4 {
    double e, x, y, z;
};
inline double Dot(const V4 &a, const V4 &b) { return a.e * b.e - a.x * b.x - a.y * b.y - a.z * b.z; }
inline V4 Sub(const V4 &a, const V4 &b) { return {a.e - b.e, a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V4 Axpy(double k, const V4 &a, const V4 &b) { return {b.e + k * a.e, b.x + k * a.x, b.y + k * a.y, b.z + k * a.z}; }
} // namespace

double OmegaEnergyDM1(double eTau, double mTau, double eRho, double eCharged, double mRho, double mCharged, double mNeutral, int nPhi) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!(eTau > mTau) || !(mRho > mCharged + mNeutral) || !(mRho < mTau) || nPhi < 4)
        return nan;
    // tau rest frame, z along the tau flight direction; the rho in the x-z plane
    // (its azimuth around z is irrelevant for S_z).
    const double gTau = eTau / mTau;
    const double bTau = std::sqrt(1. - 1. / (gTau * gTau));
    const double eRhoStar = (mTau * mTau + mRho * mRho) / (2. * mTau);
    const double pRhoStar = (mTau * mTau - mRho * mRho) / (2. * mTau);
    // E_rho(lab) = gamma (E* + beta p* cos theta)
    const double cth = std::clamp((eRho / gTau - eRhoStar) / (bTau * pRhoStar), -1., 1.);
    const double sth = std::sqrt(1. - cth * cth);
    const double g1 = eRhoStar / mRho, b1 = pRhoStar / eRhoStar;   // rho boost in the tau frame
    const double eChStar = (mRho * mRho + mCharged * mCharged - mNeutral * mNeutral) / (2. * mRho);
    const double qStar = std::sqrt(std::max(0., eChStar * eChStar - mCharged * mCharged));
    // E_charged(lab) = A + B cos(psi) - D sin(psi) cos(phi)
    const double A = gTau * g1 * eChStar * (1. + bTau * b1 * cth);
    const double B = gTau * g1 * qStar * (b1 + bTau * cth);
    const double D = gTau * bTau * qStar * sth;
    const double R = eCharged - A;

    const V4 rho{eRhoStar, pRhoStar * sth, 0., pRhoStar * cth};
    const V4 P{mTau, 0., 0., 0.};
    const V4 N = Sub(P, rho);
    const double rho2 = Dot(rho, rho);
    double sumS = 0., sumU = 0.;
    for (int k = 0; k < nPhi; ++k) {
        const double phi = (k + 0.5) * 2. * M_PI / nPhi;
        const double K = D * std::cos(phi);
        const double disc = B * B + K * K - R * R;
        if (!(disc > 0.))
            continue;
        for (const int sgn : {+1, -1}) {
            const double c = (B * R + sgn * K * std::sqrt(disc)) / (B * B + K * K);
            if (std::abs(c) > 1.)
                continue;
            const double s = std::sqrt(1. - c * c);
            if (std::abs(B * c - K * s - R) > 1e-6 * std::abs(eCharged) + 1e-9)
                continue; 
            const double jac = std::abs(B + K * c / std::max(s, 1e-12));
            // charged pion in the tau frame: rho-frame momentum (q* s cos phi, q* s sin phi, q* c)
            // in the basis (x', y, n_rho) with x' = (cos th, 0, -sin th), then boosted along n_rho.
            const double pt1 = qStar * s * std::cos(phi), pt2 = qStar * s * std::sin(phi);
            const double pPar = g1 * (qStar * c + b1 * eChStar);
            const V4 ch{g1 * (eChStar + b1 * qStar * c), pt1 * cth + pPar * sth, pt2, -pt1 * sth + pPar * cth};
            const V4 pi0 = Sub(rho, ch);
            const V4 q = Sub(ch, pi0);
            const V4 J = Axpy(-Dot(q, rho) / rho2, rho, q);
            const double JN = Dot(J, N), JP = Dot(J, P), JJ = Dot(J, J), NP = Dot(N, P);
            const double U = 2. * JN * JP - JJ * NP;
            const double Sz = mTau * (2. * JN * J.z - JJ * N.z);
            sumS += Sz / jac;
            sumU += U / jac;
        }
    }
    if (!(sumU > 0.))
        return nan;
    return std::clamp(sumS / sumU, -1., 1.);
}

double RawTemplateWeight(float omega, float pTarget, float pGen) {
    if (pTarget == pGen)
        return 1.; 
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
