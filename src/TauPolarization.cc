#include <TauPolarization/TauPolarization.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include <TauPolarization/GenTauBuilder.h>
#include <TauPolarization/TauPolarimeter.h>

using namespace TauPol;

namespace {

using TauWP = TauView::WorkingPoint;
using TauWPvsMu = TauView::WorkingPointVsMu;

const TString kMuonTrigger = "HLT_IsoMu24";
constexpr float kMuonPtMin = 26.f;
constexpr float kMuonEtaMax = 2.4f;
constexpr float kMuonIsoMax = 0.15f;
constexpr float kVetoLeptonPtMin = 10.f;
constexpr float kTauPtMin = 30.f;
constexpr float kTauEtaMax = 2.3f;
constexpr float kDRMuTauMin = 0.5f;
constexpr float kJetPtMin = 30.f;
constexpr float kJetEtaMax = 2.5f;
constexpr float kJetCleanDR = 0.4f;
constexpr int kNJetMin = 2;
constexpr int kNBJetMin = 1;
constexpr float kGenMatchDR = 0.2f;

// field order is {vsJet, vsE, vsMu, requireNewDM, maxDz, requireCalibratedDecayMode}.
// calibrated-DM requirement keeps DM 0, 1, 10, 11 and is needed for the energy scale.
const TauView::TauID kTauID{TauWP::MEDIUM, TauWP::VVLOOSE, TauWPvsMu::TIGHT, true, 0.2f, true};

const TString kMuonIDKey = "NUM_TightID_DEN_TrackerMuons";
const TString kMuonIsoKey = "NUM_TightPFIso_DEN_TightID";
const TString kMuonTriggerKey = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";

enum CutStep {
    kAll = 0, kGolden, kTrigger, kFilters, kOneMuon, kLeptonVeto, kTrigMatch,
    kTau, kJets, kBJets, kOS, kNSteps
};

// cutflow bins : 
// all, golden lumi, HLT_IsoMu24, MET filters + jet veto map,
// 1 tight iso mu (no 2nd loose mu), no loose e, mu trigger match, tau_h, >= 2 jets,
// >= 1 b jet, OS (full weight).

float TransverseMass(const TLorentzVector &l, const Particle &met) {
    const float dphi = l.DeltaPhi(met);
    return std::sqrt(std::max(0.f, 2.f * static_cast<float>(l.Pt() * met.Pt()) * (1.f - std::cos(dphi))));
}

const char *FlavName(int flav) {
    switch (flav) {
    case 5: return "GenuineTau";
    case 1: case 3: return "EtoTau";
    case 2: case 4: return "MuToTau";
    case 0: return "JetToTau";
    default: return "Data";
    }
}

} // namespace

void TauPolarization::initializeAnalyzer() {
    if (!myCorr)
        myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
    if (!HasFlag("NoTES"))
        EnableTauEnergyScale(kTauID);
    RegisterTasks();
    InitializeTasks("TauPolAll");
}

void TauPolarization::RegisterTasks() {
    RegisterTask("GenStudy", TaskOptions{true, TaskSystematicPolicy::CentralOnly},
                 [this] { ValidateGenStudy(); }, [this] { BookGenStudy(); },
                 [this] { RunGenStudy(); });
    RegisterTask("MuTauSR", TaskOptions{false, TaskSystematicPolicy::CentralOnly},
                 [] {}, [this] { BookMuTauSR(); }, [this] { RunMuTauSR(); });
}

void TauPolarization::executeEvent() {
    if (!HasTasksForSystematic(true))
        return;
    RunTasks(true);
}

// ---------------------------------------------------------------- GenStudy

void TauPolarization::ValidateGenStudy() {
    if (IsDATA)
        throw std::runtime_error("TauPolarization: GenStudy runs on MC only");
}

void TauPolarization::BookGenStudy() {
    genTaus_ = Output().Book("GenTaus");
    genTaus_.Field("weight", genRow_.weight)
        .Field("charge", genRow_.charge)
        .Field("origin", genRow_.origin)
        .Field("ancestorPdgId", genRow_.ancestorPdgId)
        .Field("genDM", genRow_.genDM)
        .Field("visTauDM", genRow_.visTauDM)
        .Field("chargedPdgId", genRow_.chargedPdgId)
        .Field("nPhotons", genRow_.nPhotons)
        .Field("valid", genRow_.valid)
        .Field("clamped", genRow_.clamped)
        .Field("omega", genRow_.omega)
        .Field("cosThetaStar", genRow_.cosThetaStar)
        .Field("cosPsi", genRow_.cosPsi)
        .Field("x", genRow_.x)
        .Field("upsilon", genRow_.upsilon)
        .Field("tauPt", genRow_.tauPt)
        .Field("tauEta", genRow_.tauEta)
        .Field("tauE", genRow_.tauE)
        .Field("visPt", genRow_.visPt)
        .Field("visEta", genRow_.visEta)
        .Field("visE", genRow_.visE)
        .Field("visM", genRow_.visM)
        .Field("visTauDR", genRow_.visTauDR)
        .Field("visTauRelDE", genRow_.visTauRelDE)
        .Field("nuResidual", genRow_.nuResidual);
}

namespace {

const char *DMName(int dm) {
    switch (dm) {
    case kDMLeptonic: return "Lep";
    case kDM0: return "DM0";
    case kDM1: return "DM1";
    case kDM2: return "DM2";
    case kDM10: return "DM10";
    case kDM11: return "DM11";
    default: return "DMOther";
    }
}

} // namespace

void TauPolarization::RunGenStudy() {
    const float weight = MCweight();
    auto hists = Hists("GenStudy");
    hists.Fill("Events", 0.5f, weight, 1, 0.f, 1.f);

    const GenViewCollection particles = GetAllGenViews();
    GenRecord rec;
    const std::size_t n = particles.size();
    rec.pt.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto p = particles[i];
        rec.pt.push_back(p.Pt());
        rec.eta.push_back(p.Eta());
        rec.phi.push_back(p.Phi());
        rec.mass.push_back(p.Mass());
        rec.pdgId.push_back(p.PdgId());
        rec.mother.push_back(p.MotherIndex());
        rec.statusFlags.push_back(p.StatusFlags());
    }
    int nOther = 0;
    const std::vector<GenTau> taus = BuildGenTausFromW(rec, &nOther);
    hists.Fill("NTausFromW", static_cast<float>(taus.size()), weight, 5, 0.f, 5.f);
    hists.Fill("NTausOtherOrigin", static_cast<float>(nOther), weight, 5, 0.f, 5.f);
    if (taus.empty())
        return;

    const GenVisTauViewCollection visTaus = GetAllGenVisTauViews();
    for (const GenTau &t : taus) {
        const PolarimeterResult r = ComputePolarimeter(t.tau, t.charged, t.neutral, t.genDM);

        GenTauRow &row = genRow_;
        row = GenTauRow{};
        row.weight = weight;
        row.charge = t.charge;
        row.origin = static_cast<int>(t.origin);
        row.ancestorPdgId = t.ancestorPdgId;
        row.genDM = t.genDM;
        row.chargedPdgId = t.chargedPdgId;
        row.nPhotons = t.nPhotons;
        row.valid = r.valid;
        row.clamped = r.clamped;
        row.omega = r.omega;
        row.cosThetaStar = r.cosThetaStar;
        row.cosPsi = r.cosPsi;
        row.x = r.x;
        row.upsilon = r.upsilon;
        row.tauPt = t.tau.Pt();
        row.tauEta = t.tau.Eta();
        row.tauE = t.tau.E();
        if (t.visible.E() > 0.) {
            row.visPt = t.visible.Pt();
            row.visEta = t.visible.Pt() > 0. ? t.visible.Eta() : 0.f;
            row.visE = t.visible.E();
            row.visM = t.visible.M();
        }
        for (std::size_t v = 0; v < visTaus.size(); ++v) {
            const auto vt = visTaus[v];
            if (vt.GenPartIdxMother() != static_cast<short>(t.index))
                continue;
            row.visTauDM = vt.Status();
            const TLorentzVector g = vt.P4();
            row.visTauDR = t.visible.Pt() > 0. ? t.visible.DeltaR(g) : -1.f;
            row.visTauRelDE = g.E() > 0. ? (t.visible.E() - g.E()) / g.E() : -9.f;
            break;
        }
        if (t.hasNu && !t.leptonic) {
            const TLorentzVector d = (t.tau - t.visible) - t.nu;
            row.nuResidual = std::sqrt(d.E() * d.E() + d.Vect().Mag2()) / t.tau.E();
        }
        genTaus_.Fill();

        // histograms  : {origin} x {charge} x {DM}.
        const std::string dir = std::string(t.origin == TauOrigin::TopW ? "Top" : "W") + "/" +
                                (t.charge < 0 ? "TauMinus" : "TauPlus") + "/" + DMName(t.genDM) + "/";
        hists.Fill(dir + "X", r.x, weight, 50, 0.f, 1.f);
        if (t.genDM == kDMLeptonic)
            continue;
        hists.Fill(dir + "Upsilon", r.upsilon, weight, 50, -1.f, 1.f);
        if (!r.valid) {
            hists.Fill(dir + "NInvalid", 0.5f, weight, 1, 0.f, 1.f);
            continue;
        }
        hists.Fill(dir + "Omega", r.omega, weight, 50, -1.f, 1.f);
        hists.Fill(dir + "CosThetaStar", r.cosThetaStar, weight, 50, -1.f, 1.f);
        if (t.genDM == kDM1) {
            hists.Fill(dir + "CosPsi", r.cosPsi, weight, 50, -1.f, 1.f);
            hists.Fill(dir + "X_vs_Upsilon", r.x, r.upsilon, weight, 25, 0.f, 1.f, 25, -1.f, 1.f);
        }
    }
}

// ---------------------------------------------------------------- MuTauSR

void TauPolarization::BookMuTauSR() {
    recoTaus_ = Output().Book("RecoTaus");
    recoTaus_.Field("weight", recoRow_.weight)
        .Field("oppositeSign", recoRow_.oppositeSign)
        .Field("nJet", recoRow_.nJet)
        .Field("nBJet", recoRow_.nBJet)
        .Field("muonPt", recoRow_.muonPt)
        .Field("muonEta", recoRow_.muonEta)
        .Field("tauPt", recoRow_.tauPt)
        .Field("tauMiniAODPt", recoRow_.tauMiniAODPt)
        .Field("tauEta", recoRow_.tauEta)
        .Field("tauCharge", recoRow_.tauCharge)
        .Field("tauDM", recoRow_.tauDM)
        .Field("tauDMUParT", recoRow_.tauDMUParT)
        .Field("leadTkRatio", recoRow_.leadTkRatio)
        .Field("upsilon", recoRow_.upsilon)
        .Field("upsilonTesUp", recoRow_.upsilonTesUp)
        .Field("upsilonTesDown", recoRow_.upsilonTesDown)
        .Field("tauPtTesUp", recoRow_.tauPtTesUp)
        .Field("tauPtTesDown", recoRow_.tauPtTesDown)
        .Field("mVis", recoRow_.mVis)
        .Field("mT", recoRow_.mT)
        .Field("met", recoRow_.met)
        .Field("genFlav", recoRow_.genFlav)
        .Field("genDM", recoRow_.genDM)
        .Field("genCharge", recoRow_.genCharge)
        .Field("genOrigin", recoRow_.genOrigin)
        .Field("genNPhotons", recoRow_.genNPhotons)
        .Field("genValid", recoRow_.genValid)
        .Field("genClamped", recoRow_.genClamped)
        .Field("genOmega", recoRow_.genOmega)
        .Field("genUpsilon", recoRow_.genUpsilon)
        .Field("genX", recoRow_.genX)
        .Field("genVisPt", recoRow_.genVisPt)
        .Field("genTauPt", recoRow_.genTauPt)
        .Field("genMatchDR", recoRow_.genMatchDR);
}

void TauPolarization::RunMuTauSR() {
    auto hists = Hists("MuTauSR");
    const Event ev = GetEvent();
    const float preWeight = IsDATA ? 1.f : MCweight() * ev.GetTriggerLumi(kMuonTrigger) * myCorr->GetPUWeight(ev.nTrueInt());
    
    const auto step = [&](int s, float w) {
        hists.Fill("Cutflow", static_cast<float>(s), w, kNSteps, 0.f, static_cast<float>(kNSteps));
        hists.Fill("CutflowRaw", static_cast<float>(s), 1.f, kNSteps, 0.f, static_cast<float>(kNSteps));
    };
    step(kAll, preWeight);
    if (IsDATA && !myCorr->IsGoldenLumi(RunNumber, ev.lumi()))
        return;
    step(kGolden, preWeight);
    if (!ev.PassTrigger(kMuonTrigger))
        return;
    step(kTrigger, preWeight);

    const JetViewCollection jets = GetAllJetViews();
    if (!PassNoiseFilterCommon(jets, ev, Event::MET_Type::PUPPI) || !PassJetVetoMap(jets))
        return;
    step(kFilters, preWeight);

    // --- exactly one tight isolated muon, no other loose muon -------------
    const MuonViewCollection muons = GetAllMuonViews();
    std::vector<std::size_t> tight;
    int nLooseMuons = 0;
    for (std::size_t i = 0; i < muons.size(); ++i) {
        const auto mu = muons[i];
        if (std::fabs(mu.Eta()) > kMuonEtaMax)
            continue;
        if (mu.Pt() > kVetoLeptonPtMin && mu.PassID(MuonView::MuonID::POG_LOOSE))
            ++nLooseMuons;
        if (mu.Pt() > kMuonPtMin && mu.PassID(MuonView::MuonID::POG_TIGHT) && mu.PfRelIso04() < kMuonIsoMax)
            tight.push_back(i);
    }
    if (tight.size() != 1 || nLooseMuons != 1)
        return;
    step(kOneMuon, preWeight);

    const ElectronViewCollection electrons = GetAllElectronViews();
    for (std::size_t i = 0; i < electrons.size(); ++i) {
        const auto e = electrons[i];
        if (e.Pt() > kVetoLeptonPtMin && e.PassID(ElectronView::ElectronID::POG_LOOSE))
            return;
    }
    step(kLeptonVeto, preWeight);

    const std::size_t muIdx = tight.front();
    const auto muon = muons[muIdx];
    const TrigObjViewCollection trigobjs = GetAllTrigObjViews();
    const auto matched = MatchViewsToTrigObjs(muons, trigobjs, 0.1f, 13);
    if (muIdx >= matched.size() || matched[muIdx] < 0)
        return;
    step(kTrigMatch, preWeight);

    const TauViewCollection taus = GetAllTauViews();
    std::size_t tauIdx = taus.size();
    for (std::size_t i = 0; i < taus.size(); ++i) {
        const auto t = taus[i];
        if (t.Pt() < kTauPtMin || std::fabs(t.Eta()) > kTauEtaMax)
            continue;
        if (!t.PassID(kTauID))
            continue;
        if (t.P4().DeltaR(muon.P4()) < kDRMuTauMin)
            continue;
        if (tauIdx == taus.size() || t.Pt() > taus[tauIdx].Pt())
            tauIdx = i;
    }
    if (tauIdx == taus.size())
        return;
    const auto tau = taus[tauIdx];
    step(kTau, preWeight);

    // jets, b-tagging
    const float bWP = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
    int nJet = 0, nBJet = 0;
    for (std::size_t i = 0; i < jets.size(); ++i) {
        const auto j = jets[i];
        if (j.Pt() < kJetPtMin || std::fabs(j.Eta()) > kJetEtaMax)
            continue;
        if (!myCorr->PassJetID(j, JetView::JetID::TIGHT))
            continue;
        if (j.P4().DeltaR(muon.P4()) < kJetCleanDR || j.P4().DeltaR(tau.P4()) < kJetCleanDR)
            continue;
        ++nJet;
        if (j.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B) > bWP)
            ++nBJet;
    }
    if (nJet < kNJetMin)
        return;
    step(kJets, preWeight);
    if (nBJet < kNBJetMin)
        return;
    step(kBJets, preWeight);

    // full weighting
    float weight = preWeight;
    if (!IsDATA) {
        weight *= myCorr->GetTauIDSF_vsJet(kTauID, tau);
        weight *= myCorr->GetTauIDSF_vsE(kTauID, tau);
        weight *= myCorr->GetTauIDSF_vsMu(kTauID, tau);
        weight *= myCorr->GetMuonIDSF(kMuonIDKey, muon);
        weight *= myCorr->GetMuonIDSF(kMuonIsoKey, muon);
        const float effData = myCorr->GetMuonTriggerEff(kMuonTriggerKey, muon.Eta(), muon.Pt(), true);
        const float effMC = myCorr->GetMuonTriggerEff(kMuonTriggerKey, muon.Eta(), muon.Pt(), false);
        if (effMC > 0.f)
            weight *= effData / effMC;
        // TODO : add in b-tagging SF
    }

    const bool os = muon.Charge() * tau.Charge() < 0;
    if (os)
        step(kOS, weight);

    Particle met = ev.GetMETVector(Event::MET_Type::PUPPI);
    PropagateTauSystToMET(taus, met);

    RecoTauRow &row = recoRow_;
    row = RecoTauRow{};
    row.weight = weight;
    row.oppositeSign = os;
    row.nJet = nJet;
    row.nBJet = nBJet;
    row.muonPt = muon.Pt();
    row.muonEta = muon.Eta();
    row.tauPt = tau.Pt();
    row.tauMiniAODPt = tau.MiniAODPt();
    row.tauEta = tau.Eta();
    row.tauCharge = tau.Charge();
    row.tauDM = static_cast<int>(tau.DecayMode());
    row.tauDMUParT = static_cast<int>(Tau_decayModeUParT[tauIdx]);
    row.leadTkRatio = Tau_leadTkPtOverTauPt[tauIdx];
    row.upsilon = UpsilonFromRatio(ShiftedRatio(row.leadTkRatio, tau.MiniAODPt(), tau.Pt()));
    row.tauPtTesUp = tau.TesPtUp();
    row.tauPtTesDown = tau.TesPtDown();
    row.upsilonTesUp = UpsilonFromRatio(ShiftedRatio(row.leadTkRatio, tau.MiniAODPt(), tau.TesPtUp()));
    row.upsilonTesDown = UpsilonFromRatio(ShiftedRatio(row.leadTkRatio, tau.MiniAODPt(), tau.TesPtDown()));
    row.mVis = (muon.P4() + tau.P4()).M();
    row.mT = TransverseMass(muon.P4(), met);
    row.met = met.Pt();
    row.genFlav = IsDATA ? -1 : static_cast<int>(tau.GenPartFlav());

    if (!IsDATA && row.genFlav == 5) {
        const GenViewCollection particles = GetAllGenViews();
        GenRecord rec;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const auto p = particles[i];
            rec.pt.push_back(p.Pt());
            rec.eta.push_back(p.Eta());
            rec.phi.push_back(p.Phi());
            rec.mass.push_back(p.Mass());
            rec.pdgId.push_back(p.PdgId());
            rec.mother.push_back(p.MotherIndex());
            rec.statusFlags.push_back(p.StatusFlags());
        }
        const std::vector<GenTau> genTaus = BuildGenTausFromW(rec, nullptr);
        const GenTau *best = nullptr;
        float bestDR = kGenMatchDR;
        for (const GenTau &g : genTaus) {
            if (g.leptonic || g.visible.Pt() <= 0.)
                continue;
            const float dr = g.visible.DeltaR(tau.P4());
            if (dr < bestDR) {
                bestDR = dr;
                best = &g;
            }
        }
        if (best) {
            const PolarimeterResult r = ComputePolarimeter(best->tau, best->charged, best->neutral, best->genDM);
            row.genDM = best->genDM;
            row.genCharge = best->charge;
            row.genOrigin = static_cast<int>(best->origin);
            row.genNPhotons = best->nPhotons;
            row.genValid = r.valid;
            row.genClamped = r.clamped;
            row.genOmega = r.omega;
            row.genUpsilon = r.upsilon;
            row.genX = r.x;
            row.genVisPt = best->visible.Pt();
            row.genTauPt = best->tau.Pt();
            row.genMatchDR = bestDR;
        }
    }
    recoTaus_.Fill();

    // histogramss  : {OS,SS} x {DM} x {flavour}.
    const std::string dir = std::string(os ? "OS" : "SS") + "/DM" + std::to_string(row.tauDM) + "/" +
                            FlavName(row.genFlav) + "/";
    hists.Fill(dir + "Upsilon", row.upsilon, weight, 40, -1.f, 1.f);
    hists.Fill(dir + "TauPt", row.tauPt, weight, 40, 0.f, 200.f);
    hists.Fill(dir + "MVis", row.mVis, weight, 40, 0.f, 200.f);
}
