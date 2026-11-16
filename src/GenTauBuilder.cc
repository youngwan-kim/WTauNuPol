#include <TauPolarization/GenTauBuilder.h>

#include <cstdlib>

namespace TauPol {

TLorentzVector GenRecord::P4(std::size_t i) const {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt[i], eta[i], phi[i], PdgMass(pdgId[i], mass[i]));
    return v;
}

namespace {

bool IsNeutrino(int a) { return a == 12 || a == 14 || a == 16; }
bool IsChargedHadron(int a) { return a == 211 || a == 321; }

int FirstDifferentAncestor(const GenRecord &gen, int i) {
    const int a = std::abs(gen.pdgId[i]);
    int k = gen.mother[i];
    while (k >= 0 && std::abs(gen.pdgId[k]) == a)
        k = gen.mother[k];
    return k;
}

} // namespace

std::vector<GenTau> BuildGenTausFromW(const GenRecord &gen, int *nOtherOrigin) {
    std::vector<GenTau> out;
    const int n = static_cast<int>(gen.size());
    for (int i = 0; i < n; ++i) {
        if (std::abs(gen.pdgId[i]) != 15 || !(gen.statusFlags[i] & kIsLastCopy))
            continue;
        const int anc = FirstDifferentAncestor(gen, i);
        if (anc < 0 || std::abs(gen.pdgId[anc]) != 24) {
            if (nOtherOrigin)
                ++*nOtherOrigin;
            continue;
        }
        GenTau t;
        t.index = i;
        t.charge = gen.pdgId[i] > 0 ? -1 : +1;
        t.ancestorPdgId = gen.pdgId[anc];
        const int wAnc = FirstDifferentAncestor(gen, anc);
        t.origin = (wAnc >= 0 && std::abs(gen.pdgId[wAnc]) == 6) ? TauOrigin::TopW : TauOrigin::W;
        t.tau = gen.P4(i);

        double leadChargedE = -1.;
        for (int j = 0; j < n; ++j) {
            if (gen.mother[j] != i)
                continue;
            const int a = std::abs(gen.pdgId[j]);
            const TLorentzVector p = gen.P4(j);
            if (IsNeutrino(a)) {
                if (a == 16) {
                    t.hasNu = true;
                    t.nu = p;
                }
                continue;
            }
            t.visible += p;
            if (a == 11 || a == 13) {
                t.leptonic = true;
            } else if (IsChargedHadron(a)) {
                ++t.nChargedHadrons;
                if (p.E() > leadChargedE) {
                    leadChargedE = p.E();
                    t.charged = p;
                    t.chargedPdgId = gen.pdgId[j];
                }
            } else if (a == 111) {
                ++t.nPi0;
            } else if (a == 22) {
                ++t.nPhotons;
            } else {
                ++t.nOtherNeutral; 
            }
        }
        if (t.leptonic) {
            t.genDM = kDMLeptonic;
        } else if (t.nOtherNeutral > 0) {
            t.genDM = kDMOther;
        } else {
            t.genDM = GenDecayModeFromCounts(t.nChargedHadrons, t.nPi0);
        }
        t.neutral = t.visible - t.charged;
        out.push_back(t);
    }
    return out;
}

} // namespace TauPol
