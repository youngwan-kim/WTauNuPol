#include <TauPolarization/GenTauBuilder.h>

#include <gtest/gtest.h>

using namespace TauPol;

namespace {

struct RecordMaker {
    GenRecord rec;
    int Add(int pdg, int mother, float pt, float eta, float phi, bool lastCopy = true, float mass = 0.f) {
        rec.pt.push_back(pt);
        rec.eta.push_back(eta);
        rec.phi.push_back(phi);
        rec.mass.push_back(mass);
        rec.pdgId.push_back(pdg);
        rec.mother.push_back(mother);
        rec.statusFlags.push_back(lastCopy ? kIsLastCopy : 0);
        return static_cast<int>(rec.pt.size()) - 1;
    }
};

} // namespace

TEST(GenTauBuilder, WalksCopiesAndClassifies) {
    RecordMaker m;
    // W- -> tau- (two copies) -> pi- pi0 nu_tau. NanoAOD stores mass 0 for all of these.
    const int w = m.Add(-24, -1, 30.f, 0.f, 0.f, true, 80.4f);
    const int tauFirst = m.Add(15, w, 40.f, 0.5f, 1.f, false);
    const int tau = m.Add(15, tauFirst, 40.f, 0.5f, 1.f);
    m.Add(-211, tau, 20.f, 0.5f, 1.f);
    m.Add(111, tau, 12.f, 0.52f, 1.02f);
    m.Add(16, tau, 8.f, 0.48f, 0.98f);
    // t -> W+ -> tau+ -> pi+ nubar
    const int top = m.Add(6, -1, 100.f, 0.f, 2.f, true, 172.5f);
    const int w2 = m.Add(24, top, 60.f, 0.1f, 2.f, true, 80.4f);
    const int tau2 = m.Add(-15, w2, 50.f, 0.1f, 2.f);
    m.Add(211, tau2, 30.f, 0.1f, 2.f);
    m.Add(-16, tau2, 20.f, 0.1f, 2.f);
    // Z -> tau tau (other origin), and a leptonic W tau.
    const int z = m.Add(23, -1, 10.f, 0.f, 0.f, true, 91.f);
    m.Add(15, z, 30.f, 1.f, 0.f);
    const int w3 = m.Add(24, -1, 20.f, -1.f, 0.f, true, 80.4f);
    const int tau3 = m.Add(-15, w3, 25.f, -1.f, 0.f);
    m.Add(-13, tau3, 10.f, -1.f, 0.f);
    m.Add(14, tau3, 5.f, -1.f, 0.f);
    m.Add(-16, tau3, 10.f, -1.f, 0.f);

    int nOther = 0;
    const auto taus = BuildGenTausFromW(m.rec, &nOther);
    ASSERT_EQ(taus.size(), 3u);
    EXPECT_EQ(nOther, 1);

    const GenTau &a = taus[0];
    EXPECT_EQ(a.index, tau);
    EXPECT_EQ(a.charge, -1);
    EXPECT_EQ(a.origin, TauOrigin::W);
    EXPECT_EQ(a.ancestorPdgId, -24);
    EXPECT_EQ(a.genDM, kDM1);
    EXPECT_EQ(a.chargedPdgId, -211);
    EXPECT_NEAR(a.tau.M(), kTauMass, 1e-4);
    EXPECT_NEAR(a.charged.M(), kPiChargedMass, 1e-4);
    EXPECT_NEAR(a.neutral.M(), kPiNeutralMass, 1e-3);
    EXPECT_TRUE(a.hasNu);

    const GenTau &b = taus[1];
    EXPECT_EQ(b.charge, +1);
    EXPECT_EQ(b.origin, TauOrigin::TopW);
    EXPECT_EQ(b.genDM, kDM0);
    EXPECT_NEAR(b.neutral.E(), 0., 1e-9);

    const GenTau &c = taus[2];
    EXPECT_TRUE(c.leptonic);
    EXPECT_EQ(c.genDM, kDMLeptonic);
}

TEST(GenTauBuilder, OtherNeutralHadronIsOtherMode) {
    RecordMaker m;
    const int w = m.Add(24, -1, 30.f, 0.f, 0.f, true, 80.4f);
    const int tau = m.Add(-15, w, 40.f, 0.f, 0.f);
    m.Add(211, tau, 20.f, 0.f, 0.f);
    m.Add(311, tau, 10.f, 0.f, 0.f, true, 0.498f);
    m.Add(-16, tau, 10.f, 0.f, 0.f);
    const auto taus = BuildGenTausFromW(m.rec);
    ASSERT_EQ(taus.size(), 1u);
    EXPECT_EQ(taus[0].genDM, kDMOther);
}

TEST(GenTauBuilder, PhotonDaughterIsOtherMode) {
    // tau -> pi pi0 gamma nu (omega -> pi0 gamma with the omega pruned) has no rho polarimeter.
    RecordMaker m;
    const int w = m.Add(24, -1, 30.f, 0.f, 0.f, true, 80.4f);
    const int tau = m.Add(-15, w, 40.f, 0.f, 0.f);
    m.Add(211, tau, 15.f, 0.f, 0.f);
    m.Add(111, tau, 10.f, 0.f, 0.f);
    m.Add(22, tau, 5.f, 0.f, 0.f);
    m.Add(-16, tau, 10.f, 0.f, 0.f);
    const auto taus = BuildGenTausFromW(m.rec);
    ASSERT_EQ(taus.size(), 1u);
    EXPECT_EQ(taus[0].nPhotons, 1);
    EXPECT_EQ(taus[0].genDM, kDMOther);
}

TEST(GenTauBuilder, AcceptZAddsZTaus) {
    // Z -> tau- tau+ (rho, pi) plus a W tau: the W-only builder sees one tau, acceptZ three.
    RecordMaker m;
    const int z = m.Add(23, -1, 10.f, 0.f, 0.f, true, 91.f);
    const int taum = m.Add(15, z, 35.f, 0.3f, 0.5f);
    m.Add(-211, taum, 18.f, 0.3f, 0.5f);
    m.Add(111, taum, 10.f, 0.32f, 0.52f);
    m.Add(16, taum, 7.f, 0.28f, 0.48f);
    const int taup = m.Add(-15, z, 30.f, -0.4f, -2.5f);
    m.Add(211, taup, 20.f, -0.4f, -2.5f);
    m.Add(-16, taup, 10.f, -0.4f, -2.5f);
    const int w = m.Add(24, -1, 30.f, 0.f, 1.f, true, 80.4f);
    const int tauw = m.Add(-15, w, 40.f, 0.f, 1.f);
    m.Add(211, tauw, 25.f, 0.f, 1.f);
    m.Add(-16, tauw, 15.f, 0.f, 1.f);

    int nOtherW = 0, nOtherZ = 0;
    const auto wOnly = BuildGenTausFromW(m.rec, &nOtherW);
    ASSERT_EQ(wOnly.size(), 1u);
    EXPECT_EQ(nOtherW, 2);
    const auto withZ = BuildGenTaus(m.rec, true, &nOtherZ);
    ASSERT_EQ(withZ.size(), 3u);
    EXPECT_EQ(nOtherZ, 0);
    EXPECT_EQ(withZ[0].origin, TauOrigin::Z);
    EXPECT_EQ(withZ[0].ancestorPdgId, 23);
    EXPECT_EQ(withZ[0].charge, -1);
    EXPECT_EQ(withZ[0].genDM, kDM1);
    EXPECT_EQ(withZ[1].origin, TauOrigin::Z);
    EXPECT_EQ(withZ[1].genDM, kDM0);
    EXPECT_EQ(withZ[2].origin, TauOrigin::W);
}
