#include <TauPolarization/TauPolarimeter.h>

#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "TRandom3.h"

using namespace TauPol;

namespace {

constexpr double kRhoMass = 0.775;
constexpr double kRhoWidth = 0.149;

TVector3 RandomDirection(TRandom3 &rng) {
    const double c = rng.Uniform(-1., 1.), phi = rng.Uniform(0., 2. * M_PI);
    const double s = std::sqrt(1. - c * c);
    return {s * std::cos(phi), s * std::sin(phi), c};
}

// momentum of each daughter in a two-body decay M -> m1 m2.
double TwoBodyMomentum(double M, double m1, double m2) {
    const double a = (M * M - (m1 + m2) * (m1 + m2)) * (M * M - (m1 - m2) * (m1 - m2));
    return a > 0. ? std::sqrt(a) / (2. * M) : 0.;
}

TLorentzVector FromMomentum(const TVector3 &p, double m) {
    TLorentzVector v;
    v.SetVectM(p, m);
    return v;
}

// Ddcay products in the tau rest frame.
struct RestDecay {
    TLorentzVector charged, neutral, nu;
};

RestDecay PiDecay(const TVector3 &dir) {
    const double p = TwoBodyMomentum(kTauMass, kPiChargedMass, 0.);
    return {FromMomentum(p * dir, kPiChargedMass), TLorentzVector(), FromMomentum(-p * dir, 0.)};
}

// tau -> rho nu, rho (mass mRho) along rhoDir, pi+- along piDir in the rho frame.
RestDecay RhoDecay(double mRho, const TVector3 &rhoDir, const TVector3 &piDir) {
    const double pRho = TwoBodyMomentum(kTauMass, mRho, 0.);
    const TLorentzVector rho = FromMomentum(pRho * rhoDir, mRho);
    const double pPi = TwoBodyMomentum(mRho, kPiChargedMass, kPiNeutralMass);
    TLorentzVector c = FromMomentum(pPi * piDir, kPiChargedMass);
    TLorentzVector n = FromMomentum(-pPi * piDir, kPiNeutralMass);
    c.Boost(rho.BoostVector());
    n.Boost(rho.BoostVector());
    return {c, n, FromMomentum(-pRho * rhoDir, 0.)};
}

// put the tau at lab momentum pLab along labDir, returns lab four-vectors.
struct LabDecay {
    TLorentzVector tau, charged, neutral;
};

LabDecay ToLab(const RestDecay &rest, const TVector3 &labDir, double pLab) {
    TLorentzVector tau = FromMomentum(pLab * labDir, kTauMass);
    LabDecay out{tau, rest.charged, rest.neutral};
    out.charged.Boost(tau.BoostVector());
    if (rest.neutral.E() > 0.)
        out.neutral.Boost(tau.BoostVector());
    return out;
}

double RandomRhoMass(TRandom3 &rng) {
    while (true) {
        const double m = rng.BreitWigner(kRhoMass, kRhoWidth);
        if (m > kPiChargedMass + kPiNeutralMass + 1e-3 && m < kTauMass - 1e-3)
            return m;
    }
}

} // namespace

// 1. DM0: omega = cos(theta*) of the pion in the helicity frame.
TEST(TauPolarimeter, DM0OmegaIsCosThetaStar) {
    TRandom3 rng(1);
    for (int i = 0; i < 10000; ++i) {
        const TVector3 dir = RandomDirection(rng), labDir = RandomDirection(rng);
        const LabDecay lab = ToLab(PiDecay(dir), labDir, rng.Uniform(5., 500.));
        const PolarimeterResult r = ComputePolarimeter(lab.tau, lab.charged, lab.neutral, kDM0);
        ASSERT_TRUE(r.valid);
        TLorentzVector pi = lab.charged;
        pi.Boost(-lab.tau.BoostVector());
        const double expected = pi.Vect().Unit().Dot(lab.tau.Vect().Unit());
        EXPECT_NEAR(r.omega, expected, 1e-6);
        EXPECT_NEAR(r.omega, dir.Dot(labDir), 1e-5);
        EXPECT_NEAR(r.cosThetaStar, expected, 1e-6);
        EXPECT_NEAR(r.x, lab.charged.E() / lab.tau.E(), 1e-6);
        EXPECT_NEAR(r.upsilon, 1.f, 1e-6);
    }
}

// DM0 omega must not depend on the lab opening angle, but NanoAOD's reduced GenPart precision smears by ~1e-3. 
// Rotating the pion by 2 mrad at fixed energy leaves omega unchanged, exact inputs are never clamped.
TEST(TauPolarimeter, DM0RobustToAngularSmearing) {
    TRandom3 rng(11);
    for (int i = 0; i < 2000; ++i) {
        const TVector3 dir = RandomDirection(rng), labDir = RandomDirection(rng);
        const LabDecay lab = ToLab(PiDecay(dir), labDir, rng.Uniform(100., 1500.));
        const auto exact = ComputePolarimeter(lab.tau, lab.charged, lab.neutral, kDM0);
        ASSERT_TRUE(exact.valid);
        EXPECT_FALSE(exact.clamped);
        TLorentzVector smeared = lab.charged;
        smeared.Rotate(2e-3, lab.charged.Vect().Orthogonal());
        const auto r = ComputePolarimeter(lab.tau, smeared, lab.neutral, kDM0);
        ASSERT_TRUE(r.valid);
        EXPECT_NEAR(r.omega, exact.omega, 1e-5);
    }
}

// 2. DM1: |omega| <= 1, omega independent of the tau lab energy.
TEST(TauPolarimeter, DM1BoundedAndBoostInvariant) {
    TRandom3 rng(2);
    for (int i = 0; i < 100000; ++i) {
        const RestDecay rest = RhoDecay(RandomRhoMass(rng), RandomDirection(rng), RandomDirection(rng));
        const TVector3 labDir = RandomDirection(rng);
        const LabDecay a = ToLab(rest, labDir, rng.Uniform(1., 30.));
        const PolarimeterResult ra = ComputePolarimeter(a.tau, a.charged, a.neutral, kDM1);
        ASSERT_TRUE(ra.valid);
        ASSERT_LE(std::abs(ra.omega), 1.f);
        if (i % 10 == 0) {
            const LabDecay b = ToLab(rest, labDir, rng.Uniform(100., 1000.));
            const PolarimeterResult rb = ComputePolarimeter(b.tau, b.charged, b.neutral, kDM1);
            ASSERT_TRUE(rb.valid);
            EXPECT_NEAR(ra.omega, rb.omega, 1e-4);
            EXPECT_NEAR(ra.cosPsi, rb.cosPsi, 1e-4);
        }
    }
}

// 3. Reflecting all decay momenta in the tau rest frame flips omega.
TEST(TauPolarimeter, ReflectionFlipsOmega) {
    TRandom3 rng(3);
    auto reflect = [](const TLorentzVector &v) {
        return TLorentzVector(-v.Px(), -v.Py(), -v.Pz(), v.E());
    };
    for (int i = 0; i < 10000; ++i) {
        const RestDecay rest = RhoDecay(RandomRhoMass(rng), RandomDirection(rng), RandomDirection(rng));
        const RestDecay mirrored{reflect(rest.charged), reflect(rest.neutral), reflect(rest.nu)};
        const TVector3 labDir = RandomDirection(rng);
        const double p = rng.Uniform(5., 200.);
        const LabDecay a = ToLab(rest, labDir, p), b = ToLab(mirrored, labDir, p);
        const auto ra = ComputePolarimeter(a.tau, a.charged, a.neutral, kDM1);
        const auto rb = ComputePolarimeter(b.tau, b.charged, b.neutral, kDM1);
        ASSERT_TRUE(ra.valid && rb.valid);
        EXPECT_NEAR(ra.omega, -rb.omega, 1e-4);
        EXPECT_NEAR(ra.u, rb.u, 1e-6 * std::abs(ra.u));
    }
}

// 4. At E_tau = 1 TeV the lab energy fraction equals the rest-frame light-cone fraction.
TEST(TauPolarimeter, LightConeLimit) {
    TRandom3 rng(4);
    for (int i = 0; i < 10000; ++i) {
        const RestDecay rest = RhoDecay(RandomRhoMass(rng), RandomDirection(rng), RandomDirection(rng));
        const TVector3 labDir = RandomDirection(rng);
        const LabDecay lab = ToLab(rest, labDir, 1000.);
        const auto r = ComputePolarimeter(lab.tau, lab.charged, lab.neutral, kDM1);
        auto plus = [&](const TLorentzVector &v) { return v.E() + v.Vect().Dot(labDir); };
        const double qPlus = plus(rest.charged) - plus(rest.neutral);
        const double QPlus = plus(rest.charged) + plus(rest.neutral);
        EXPECT_NEAR(r.upsilon, qPlus / QPlus, 1e-3);
        EXPECT_NEAR(r.x, QPlus / kTauMass, 1e-3);
    }
}

// 5. Toy tau -> rho nu with the closed form and a constant-width BW. Events are
// generated flat in angles and weighted by the full rate:
//   |BW(s)|^2 x (two-body phase spaces) x U (1 + P omega).
// Expected (theory note): alpha_rho = 0.436, f_L = 0.718, <omega> = P / 3.
TEST(TauPolarimeter, ToyRhoReferenceNumbers) {
    TRandom3 rng(5);
    const double mMin = kPiChargedMass + kPiNeutralMass, mMax = kTauMass;
    const double hw = 0.5 * kRhoWidth;
    const double aMin = std::atan((mMin - kRhoMass) / hw), aMax = std::atan((mMax - kRhoMass) / hw);
    const TVector3 labDir(0., 0., 1.);

    // accumulate weighted sums for P = +1 and P = -1 from the same sample.
    double sw[2] = {0, 0}, sw2[2] = {0, 0}, sCos[2] = {0, 0}, sCos2Psi[2] = {0, 0}, sOmega[2] = {0, 0};
    const int nEvents = 400000;
    for (int i = 0; i < nEvents; ++i) {
        const double m = kRhoMass + hw * std::tan(rng.Uniform(aMin, aMax));
        const double s = m * m;
        const double cauchy = 1. / ((m - kRhoMass) * (m - kRhoMass) + hw * hw);
        const double bw2 = 1. / ((s - kRhoMass * kRhoMass) * (s - kRhoMass * kRhoMass) + kRhoMass * kRhoMass * kRhoWidth * kRhoWidth);
        const double phaseSpace = TwoBodyMomentum(kTauMass, m, 0.) / kTauMass * TwoBodyMomentum(m, kPiChargedMass, kPiNeutralMass) / m;
        const double wMass = bw2 * 2. * m / cauchy * phaseSpace;
        const LabDecay lab = ToLab(RhoDecay(m, RandomDirection(rng), RandomDirection(rng)), labDir, 20.);
        const auto r = ComputePolarimeter(lab.tau, lab.charged, lab.neutral, kDM1);
        ASSERT_TRUE(r.valid);
        for (int k = 0; k < 2; ++k) {
            const double P = k == 0 ? 1. : -1.;
            const double w = wMass * r.u * (1. + P * r.omega);
            sw[k] += w;
            sw2[k] += w * w;
            sCos[k] += w * r.cosThetaStar;
            sCos2Psi[k] += w * r.cosPsi * r.cosPsi;
            sOmega[k] += w * r.omega;
        }
    }
    for (int k = 0; k < 2; ++k) {
        const double P = k == 0 ? 1. : -1.;
        const double alpha = 3. * (sCos[k] / sw[k]) / P;
        const double fL = (sCos2Psi[k] / sw[k] - 0.2) / 0.4;
        const double omegaOverP3 = (sOmega[k] / sw[k]) / (P / 3.);
        const double nEff = sw[k] * sw[k] / sw2[k];
        RecordProperty(k == 0 ? "alpha_Pplus" : "alpha_Pminus", std::to_string(alpha));
        RecordProperty(k == 0 ? "fL_Pplus" : "fL_Pminus", std::to_string(fL));
        RecordProperty(k == 0 ? "omegaRatio_Pplus" : "omegaRatio_Pminus", std::to_string(omegaOverP3));
        RecordProperty(k == 0 ? "nEff_Pplus" : "nEff_Pminus", std::to_string(nEff));
        std::printf("[toy rho] P=%+.0f alpha=%.4f fL=%.4f <omega>/(P/3)=%.4f nEff=%.0f\n", P, alpha, fL, omegaOverP3, nEff);
        EXPECT_NEAR(alpha, 0.436, 0.01);
        EXPECT_NEAR(fL, 0.718, 0.01);
        EXPECT_NEAR(omegaOverP3, 1., 0.02);
    }
}

TEST(TauPolarimeter, ToyPionMeanOmega) {
    TRandom3 rng(6);
    const TVector3 labDir(0., 0., 1.);
    for (const double P : {1., -1.}) {
        double sw = 0., so = 0.;
        for (int i = 0; i < 200000; ++i) {
            const LabDecay lab = ToLab(PiDecay(RandomDirection(rng)), labDir, 30.);
            const auto r = ComputePolarimeter(lab.tau, lab.charged, lab.neutral, kDM0);
            const double w = 1. + P * r.omega;
            sw += w;
            so += w * r.omega;
        }
        EXPECT_NEAR(so / sw, P / 3., 0.005);
    }
}

// 6. Template weight.
TEST(TauPolarimeter, TemplateWeight) {
    for (float omega : {-1.f, -0.5f, 0.f, 0.3f, 0.99f})
        for (float p : {-1.f, -0.9f, 0.f, 0.5f})
            EXPECT_FLOAT_EQ(TemplateWeight(omega, p, p, 50.f), 1.f);
    EXPECT_FLOAT_EQ(TemplateWeight(1.f, -1.f, -1.f, 50.f), 1.f);  // 0/0 at the boundary
    EXPECT_FLOAT_EQ(TemplateWeight(0.5f, 1.f, -1.f, 50.f), 3.f);
    EXPECT_FLOAT_EQ(TemplateWeight(0.999f, 1.f, -1.f, 50.f), 50.f);
    EXPECT_FLOAT_EQ(TemplateWeight(1.f, 1.f, -1.f, 50.f), 50.f);
    EXPECT_TRUE(std::isinf(RawTemplateWeight(1.f, 1.f, -1.f)));
    EXPECT_FLOAT_EQ(TemplateWeight(0.5f, -1.f, -1.f, 50.f), 1.f);
    EXPECT_FLOAT_EQ(TemplateWeight(0.5f, -1.f, 1.f, 50.f), 1.f / 3.f);
}

TEST(TauPolarimeter, RatioHelpers) {
    EXPECT_FLOAT_EQ(UpsilonFromRatio(0.5f), 0.f);
    EXPECT_FLOAT_EQ(UpsilonFromRatio(1.f), 1.f);
    // A +3% TES shift moves the denominator, not the track.
    EXPECT_FLOAT_EQ(ShiftedRatio(0.6f, 40.f, 41.2f), 0.6f * 40.f / 41.2f);
    EXPECT_FLOAT_EQ(ShiftedRatio(0.6f, 40.f, 0.f), 0.6f);
}

TEST(TauPolarimeter, DecayModeAndMasses) {
    EXPECT_EQ(GenDecayModeFromCounts(1, 0), 0);
    EXPECT_EQ(GenDecayModeFromCounts(1, 1), 1);
    EXPECT_EQ(GenDecayModeFromCounts(1, 2), 2);
    EXPECT_EQ(GenDecayModeFromCounts(1, 3), 15);
    EXPECT_EQ(GenDecayModeFromCounts(3, 0), 10);
    EXPECT_EQ(GenDecayModeFromCounts(3, 1), 11);
    EXPECT_EQ(GenDecayModeFromCounts(5, 0), 15);
    EXPECT_DOUBLE_EQ(PdgMass(-15, 0.), kTauMass);
    EXPECT_DOUBLE_EQ(PdgMass(211, 0.), kPiChargedMass);
    EXPECT_DOUBLE_EQ(PdgMass(111, 0.), kPiNeutralMass);
    EXPECT_DOUBLE_EQ(PdgMass(221, 0.5469), 0.5469);
}

TEST(TauPolarimeter, InvalidModes) {
    TLorentzVector tau, c;
    tau.SetPtEtaPhiM(40., 0.3, 1., kTauMass);
    c.SetPtEtaPhiM(20., 0.3, 1., kPiChargedMass);
    for (int dm : {kDMLeptonic, kDM2, kDM10, kDM11, kDMOther}) {
        const auto r = ComputePolarimeter(tau, c, TLorentzVector(), dm);
        EXPECT_FALSE(r.valid);
        EXPECT_EQ(r.genDM, dm);
        EXPECT_NEAR(r.x, c.E() / tau.E(), 1e-6);
    }
    // DM1 with no neutral system cannot form the current.
    EXPECT_FALSE(ComputePolarimeter(tau, c, TLorentzVector(), kDM1).valid);
}

// Section 9 Fisher information.
TEST(TauPolarimeter, ShapeFisherInformation) {
    // Identical shapes: no shape information.
    EXPECT_NEAR(ShapeFisherInformation({1, 2, 3}, {2, 4, 6}, 0.), 0., 1e-12);
    // Perfectly separated, equal acceptance: I(P) = 1 / (1 - P^2).
    EXPECT_NEAR(ShapeFisherInformation({1, 0}, {0, 1}, 0.), 1., 1e-12);
    EXPECT_NEAR(ShapeFisherInformation({1, 0}, {0, 1}, -0.5), 1. / 0.75, 1e-12);
    // Pion cos(theta*) with fine bins: I -> integral of x^2/(1+Px) dx/2 = 1/3 at P = 0.
    std::vector<double> plus, minus;
    const int nb = 2000;
    for (int b = 0; b < nb; ++b) {
        const double x = -1. + (b + 0.5) * 2. / nb;
        plus.push_back((1. + x) / nb);
        minus.push_back((1. - x) / nb);
    }
    EXPECT_NEAR(ShapeFisherInformation(plus, minus, 0.), 1. / 3., 1e-4);
}
