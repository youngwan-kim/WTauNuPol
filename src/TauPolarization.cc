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
constexpr float kSubMuonPtMin = 20.f;    // Z -> mu mu: subleading muon
constexpr float kZMassMin = 70.f, kZMassMax = 110.f;
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
// Anti-isolated control region (plan B4): fails vsJet MEDIUM, passes VVLOOSE.
const TauView::TauID kTauIDLoose{TauWP::VVLOOSE, TauWP::VVLOOSE, TauWPvsMu::TIGHT, true, 0.2f, true};
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
    doSyst_ = HasFlag("Syst");
    isPowhegTT_ = !IsDATA && (MCSample.BeginsWith("TTLL_powheg") || MCSample.BeginsWith("TTLJ_powheg"));
    RegisterTasks();
    InitializeTasks("TauPolAll");
}

void TauPolarization::RegisterTasks() {
    RegisterTask("GenStudy", TaskOptions{true, TaskSystematicPolicy::CentralOnly},
                 [this] { ValidateGenStudy(); }, [this] { BookGenStudy(); },
                 [this] { RunGenStudy(); });
    RegisterTask("MuTauSR", TaskOptions{false, TaskSystematicPolicy::CentralOnly},
                 [] {}, [this] { BookMuTauSR(); }, [this] { RunMuTau(MuTauMode::SR); });
    RegisterTask("ZMuMu", TaskOptions{false, TaskSystematicPolicy::CentralOnly},
                 [] {}, [this] { BookZMuMu(); }, [this] { RunZMuMu(); });
    RegisterTask("MuTauCR", TaskOptions{false, TaskSystematicPolicy::CentralOnly},
                 [] {}, [this] { BookMuTauCR(); }, [this] { RunMuTau(MuTauMode::CR); });
    RegisterTask("ZTauTau", TaskOptions{false, TaskSystematicPolicy::CentralOnly},
                 [] {}, [this] { BookZTauTau(); }, [this] { RunMuTau(MuTauMode::ZTT); });
}

void TauPolarization::executeEvent() {
    if (!HasTasksForSystematic(true))
        return;
    if (doSyst_ && !IsDATA)
        ComputeModelWeights();
    RunTasks(true);
}

// SystSumW bins: 0 nominal, 1-6 muR up / down, muF up / down, muR = muF up / down,
// 7-10 ISR up / down, FSR up / down, 11 top pT, 12-13 hdamp up / down,
// 14 + i PDF member i (LHEPdfWeight[i], up to kMaxPdf). 
void TauPolarization::ComputeModelWeights() {
    using var = MyCorrection::variation;
    ModelWeights &m = model_;
    m = ModelWeights{};
    try {
        m.muRUp = GetScaleVariation(var::nom, var::up);
        m.muRDown = GetScaleVariation(var::nom, var::down);
        m.muFUp = GetScaleVariation(var::up, var::nom);
        m.muFDown = GetScaleVariation(var::down, var::nom);
        m.muRFUp = GetScaleVariation(var::up, var::up);
        m.muRFDown = GetScaleVariation(var::down, var::down);
    } catch (const std::exception &) {
        m.muRUp = m.muRDown = m.muFUp = m.muFDown = m.muRFUp = m.muRFDown = 1.f;
    }
    try {
        m.isrUp = GetPSWeight(var::up, var::nom);
        m.isrDown = GetPSWeight(var::down, var::nom);
        m.fsrUp = GetPSWeight(var::nom, var::up);
        m.fsrDown = GetPSWeight(var::nom, var::down);
    } catch (const std::exception &) {
        m.isrUp = m.isrDown = m.fsrUp = m.fsrDown = 1.f;
    }
    if (LHEPdfWeight.valid()) {
        const auto &v = LHEPdfWeight.values();
        m.pdf.assign(v.begin(), v.end());
    }
    if (isPowhegTT_) {
        // first / last copy of t and tbar (status flags bits 12 / 13, as ModellingPatch)
        constexpr unsigned short kFirstCopy = 1u << 12, kLastCopy = 1u << 13;
        const GenViewCollection gens = GetAllGenViews();
        int ft = -1, fa = -1, lt = -1, la = -1;
        for (std::size_t i = 0; i < gens.size(); ++i) {
            const int pdg = gens[i].PdgId();
            const auto flags = gens[i].StatusFlags();
            if (pdg == 6 && (flags & kFirstCopy) && ft < 0) ft = static_cast<int>(i);
            if (pdg == -6 && (flags & kFirstCopy) && fa < 0) fa = static_cast<int>(i);
            if (pdg == 6 && (flags & kLastCopy)) lt = static_cast<int>(i);
            if (pdg == -6 && (flags & kLastCopy)) la = static_cast<int>(i);
        }
        if (ft >= 0 && fa >= 0 && lt >= 0 && la >= 0) {
            m.topPt = myCorr->GetTopPtReweight(gens[lt].P4(), gens[la].P4());
            m.hdampUp = myCorr->GethDampReweight(gens[ft].P4(), gens[fa].P4(), var::up);
            m.hdampDown = myCorr->GethDampReweight(gens[ft].P4(), gens[fa].P4(), var::down);
        }
    }
    constexpr int kFixed = 14, kMaxPdf = 110, kBins = kFixed + kMaxPdf;
    const float w = MCweight();
    const float fixed[kFixed] = {1.f, m.muRUp, m.muRDown, m.muFUp, m.muFDown, m.muRFUp, m.muRFDown,
                                 m.isrUp, m.isrDown, m.fsrUp, m.fsrDown, m.topPt, m.hdampUp, m.hdampDown};
    auto h = Hists();
    for (int i = 0; i < kFixed; ++i)
        h.Fill("SystSumW", i + 0.5f, w * fixed[i], kBins, 0.f, kBins);
    for (int i = 0; i < static_cast<int>(m.pdf.size()) && i < kMaxPdf; ++i)
        h.Fill("SystSumW", kFixed + i + 0.5f, w * m.pdf[i], kBins, 0.f, kBins);
}

//  GenStudy

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
        .Field("validE", genRow_.validE)
        .Field("omegaE", genRow_.omegaE)
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
        row.validE = r.validE;
        row.omegaE = r.omegaE;
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
    BindRecoFields(recoTaus_);
}

void TauPolarization::BookMuTauCR() {
    recoTausCR_ = Output().Book("RecoTausCR");
    BindRecoFields(recoTausCR_);
}

void TauPolarization::BookZTauTau() {
    recoTausZ_ = Output().Book("ZTauTauEvents");
    BindRecoFields(recoTausZ_);
}

void TauPolarization::BindRecoFields(RNTupleHandle &h) {
    if (doSyst_) {
        RecoTauRow &r = recoRow_;
        h.Field("nJetJESUp", r.nJetJESUp).Field("nJetJESDown", r.nJetJESDown)
            .Field("nJetJERUp", r.nJetJERUp).Field("nJetJERDown", r.nJetJERDown)
            .Field("nBJetJESUp", r.nBJetJESUp).Field("nBJetJESDown", r.nBJetJESDown)
            .Field("nBJetJERUp", r.nBJetJERUp).Field("nBJetJERDown", r.nBJetJERDown)
            .Field("wTauIDvsJetUp", r.wTauIDvsJetUp).Field("wTauIDvsJetDown", r.wTauIDvsJetDown)
            .Field("wTauIDvsEUp", r.wTauIDvsEUp).Field("wTauIDvsEDown", r.wTauIDvsEDown)
            .Field("wTauIDvsMuUp", r.wTauIDvsMuUp).Field("wTauIDvsMuDown", r.wTauIDvsMuDown)
            .Field("wMuIDUp", r.wMuIDUp).Field("wMuIDDown", r.wMuIDDown)
            .Field("wMuIsoUp", r.wMuIsoUp).Field("wMuIsoDown", r.wMuIsoDown)
            .Field("wMuTrigUp", r.wMuTrigUp).Field("wMuTrigDown", r.wMuTrigDown)
            .Field("wPUUp", r.wPUUp).Field("wPUDown", r.wPUDown)
            .Field("wBtagHeavyUp", r.wBtagHeavyUp).Field("wBtagHeavyDown", r.wBtagHeavyDown)
            .Field("wBtagLightUp", r.wBtagLightUp).Field("wBtagLightDown", r.wBtagLightDown)
            .Field("wISRUp", r.wISRUp).Field("wISRDown", r.wISRDown)
            .Field("wFSRUp", r.wFSRUp).Field("wFSRDown", r.wFSRDown)
            .Field("wMuRUp", r.wMuRUp).Field("wMuRDown", r.wMuRDown)
            .Field("wMuFUp", r.wMuFUp).Field("wMuFDown", r.wMuFDown)
            .Field("wMuRFUp", r.wMuRFUp).Field("wMuRFDown", r.wMuRFDown)
            .Field("wTopPt", r.wTopPt).Field("wHdampUp", r.wHdampUp).Field("wHdampDown", r.wHdampDown)
            .Field("wPdf", r.wPdf)
            .SetHalfPrecision("wPdf");
    }
    h.Field("weight", recoRow_.weight)
        .Field("btagWeight", recoRow_.btagWeight)
        .Field("region", recoRow_.region)
        .Field("tauPassMedium", recoRow_.tauPassMedium)
        .Field("oppositeSign", recoRow_.oppositeSign)
        .Field("nJet", recoRow_.nJet)
        .Field("nBJet", recoRow_.nBJet)
        .Field("muonPt", recoRow_.muonPt)
        .Field("muonEta", recoRow_.muonEta)
        .Field("muonPhi", recoRow_.muonPhi)
        .Field("tauPt", recoRow_.tauPt)
        .Field("tauMiniAODPt", recoRow_.tauMiniAODPt)
        .Field("tauEta", recoRow_.tauEta)
        .Field("tauPhi", recoRow_.tauPhi)
        .Field("tauMass", recoRow_.tauMass)
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
        .Field("metPhi", recoRow_.metPhi)
        .Field("jet1Pt", recoRow_.jet1Pt)
        .Field("lheHT", recoRow_.lheHT)
        .Field("genFlav", recoRow_.genFlav)
        .Field("genDM", recoRow_.genDM)
        .Field("genCharge", recoRow_.genCharge)
        .Field("genOrigin", recoRow_.genOrigin)
        .Field("genNPhotons", recoRow_.genNPhotons)
        .Field("genValid", recoRow_.genValid)
        .Field("genClamped", recoRow_.genClamped)
        .Field("genOmega", recoRow_.genOmega)
        .Field("genValidE", recoRow_.genValidE)
        .Field("genOmegaE", recoRow_.genOmegaE)
        .Field("genUpsilon", recoRow_.genUpsilon)
        .Field("genX", recoRow_.genX)
        .Field("genVisPt", recoRow_.genVisPt)
        .Field("genTauPt", recoRow_.genTauPt)
        .Field("genMatchDR", recoRow_.genMatchDR);
}

// Shared mu tau_h selection 
//   SR  (MuTauSR): medium vsJet tau, >= 2 jets, >= 1 b jet -> RecoTaus (region 0).
//   CR  (MuTauCR): leading VVLoose vsJet tau, >= 2 jets, then
//       region 1: tau fails Medium, >= 1 b jet (anti iso, fake enriched)
//       region 2: tau passes Medium, 0 b jet  (W + jets enriched)
//       region 3: tau fails Medium,  0 b jet  (fake factor fail side)    -> RecoTausCR.
//   ZTT  (ZTauTau): leading VVLoose vsJet tau, 0 b jets, no jet-count requirement (applied
//       offline), region 4: tau passes Medium, region 5: fails Medium -> ZTauTauEvents;
//       gen matching also accepts taus from a Z.
// Same-sign events are kept in all trees (oppositeSign = false).
void TauPolarization::RunMuTau(MuTauMode mode) {
    const bool cr = mode == MuTauMode::CR;
    const bool ztt = mode == MuTauMode::ZTT;
    auto hists = Hists(ztt ? "ZTauTau" : (cr ? "MuTauCR" : "MuTauSR"));
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
        if (std::max({t.Pt(), t.TesPtUp(), t.TesPtDown()}) < kTauPtMin || std::fabs(t.Eta()) > kTauEtaMax)
            continue;
        if (!t.PassID(mode != MuTauMode::SR ? kTauIDLoose : kTauID))
            continue;
        if (t.P4().DeltaR(muon.P4()) < kDRMuTauMin)
            continue;
        if (tauIdx == taus.size() || t.Pt() > taus[tauIdx].Pt())
            tauIdx = i;
    }
    if (tauIdx == taus.size())
        return;
    const auto tau = taus[tauIdx];
    const bool tauMedium = tau.PassID(kTauID);
    step(kTau, preWeight);

    // jets, b-tagging
    const float bWP = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
    int nJet = 0, nBJet = 0;
    // fixed-WP b-tag event weight (BTV method 1a) over the selected jets:
    float btagWeight = 1.f;
    using var = MyCorrection::variation;
    float bHUp = 1.f, bHDown = 1.f, bLUp = 1.f, bLDown = 1.f;
    float jet1Pt = 0.f;
    for (std::size_t i = 0; i < jets.size(); ++i) {
        const auto j = jets[i];
        if (j.Pt() < kJetPtMin || std::fabs(j.Eta()) >= kJetEtaMax) 
            continue;
        if (!myCorr->PassJetID(j, JetView::JetID::TIGHT))
            continue;
        if (j.P4().DeltaR(muon.P4()) < kJetCleanDR || j.P4().DeltaR(tau.P4()) < kJetCleanDR)
            continue;
        ++nJet;
        jet1Pt = std::max(jet1Pt, static_cast<float>(j.Pt()));
        const bool tagged = j.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B) > bWP;
        if (tagged)
            ++nBJet;
        if (!IsDATA) {
            const int flav = static_cast<int>(j.HadronFlavour());
            const float sf = myCorr->GetBTaggingSF(j.Eta(), j.Pt(), flav, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
            const float eff = tagged ? 1.f : myCorr->GetBTaggingEff(j.Eta(), j.Pt(), flav, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
            const auto factor = [&](float s) {
                if (tagged)
                    return s;
                return (eff > 0.f && eff < 1.f) ? (1.f - s * eff) / (1.f - eff) : 1.f;
            };
            btagWeight *= factor(sf);
            if (doSyst_) {
                const bool heavy = flav == 4 || flav == 5;
                const float sfUp = myCorr->GetBTaggingSF(j.Eta(), j.Pt(), flav, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium, var::up);
                const float sfDown = myCorr->GetBTaggingSF(j.Eta(), j.Pt(), flav, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium, var::down);
                bHUp *= factor(heavy ? sfUp : sf);
                bHDown *= factor(heavy ? sfDown : sf);
                bLUp *= factor(heavy ? sf : sfUp);
                bLDown *= factor(heavy ? sf : sfDown);
            }
        }
    }

    const auto regionFor = [&](int nj, int nb) -> int {
        if (ztt)
            return nb == 0 ? (tauMedium ? 4 : 5) : -1;
        if (nj < kNJetMin)
            return -1;
        if (!cr)
            return nb >= kNBJetMin ? 0 : -1;
        if (!tauMedium && nb >= kNBJetMin)
            return 1;
        if (tauMedium && nb == 0)
            return 2;
        if (!tauMedium && nb == 0)
            return 3;
        return -1;
    };
    const int region = regionFor(nJet, nBJet);

    int laneNJ[4] = {nJet, nJet, nJet, nJet}, laneNB[4] = {nBJet, nBJet, nBJet, nBJet};
    bool anyLane = region >= 0;
    if (doSyst_ && !IsDATA) {
        const var jesV[4] = {var::up, var::down, var::nom, var::nom};
        const var jerV[4] = {var::nom, var::nom, var::up, var::down};
        for (int l = 0; l < 4; ++l) {
            laneNJ[l] = laneNB[l] = 0;
            for (const std::size_t i : SelectJetIndices(jets, JetView::JetID::TIGHT, kJetPtMin, kJetEtaMax, jesV[l], jerV[l])) {
                const auto j = jets[i];
                if (std::fabs(j.Eta()) >= kJetEtaMax)
                    continue;
                if (j.P4().DeltaR(muon.P4()) < kJetCleanDR || j.P4().DeltaR(tau.P4()) < kJetCleanDR)
                    continue;
                ++laneNJ[l];
                if (j.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B) > bWP)
                    ++laneNB[l];
            }
            anyLane = anyLane || regionFor(laneNJ[l], laneNB[l]) >= 0;
        }
    }
    if (ztt || nJet >= kNJetMin)
        step(kJets, preWeight);
    if (!anyLane)
        return;
    if (region >= 0)
        step(kBJets, preWeight);

    // full weighting
    float weight = preWeight;
    if (!IsDATA) {
        if (tauMedium) {
            weight *= myCorr->GetTauIDSF_vsJet(kTauID, tau);
            weight *= myCorr->GetTauIDSF_vsE(kTauID, tau);
            weight *= myCorr->GetTauIDSF_vsMu(kTauID, tau);
        } else {
            weight *= myCorr->GetTauIDSF_vsE(kTauID, tau);
            weight *= myCorr->GetTauIDSF_vsMu(kTauID, tau);
        }
        weight *= myCorr->GetMuonIDSF(kMuonIDKey, muon);
        weight *= myCorr->GetMuonIDSF(kMuonIsoKey, muon);
        const float effData = myCorr->GetMuonTriggerEff(kMuonTriggerKey, muon.Eta(), muon.Pt(), true);
        const float effMC = myCorr->GetMuonTriggerEff(kMuonTriggerKey, muon.Eta(), muon.Pt(), false);
        if (effMC > 0.f)
            weight *= effData / effMC;
        weight *= btagWeight;
    }

    const bool os = muon.Charge() * tau.Charge() < 0;
    if (os && region >= 0)
        step(kOS, weight);

    Particle met = ev.GetMETVector(Event::MET_Type::PUPPI);
    PropagateTauSystToMET(taus, met);

    RecoTauRow &row = recoRow_;
    row = RecoTauRow{};
    row.weight = weight;
    row.btagWeight = btagWeight;
    row.region = region;
    row.nJetJESUp = laneNJ[0];
    row.nJetJESDown = laneNJ[1];
    row.nJetJERUp = laneNJ[2];
    row.nJetJERDown = laneNJ[3];
    row.nBJetJESUp = laneNB[0];
    row.nBJetJESDown = laneNB[1];
    row.nBJetJERUp = laneNB[2];
    row.nBJetJERDown = laneNB[3];
    if (doSyst_ && !IsDATA) {
        const auto ratio = [](float nom, float v) { return nom > 0.f ? v / nom : 1.f; };
        if (tauMedium) {
            const float n = myCorr->GetTauIDSF_vsJet(kTauID, tau);
            row.wTauIDvsJetUp = ratio(n, myCorr->GetTauIDSF_vsJet(kTauID, tau, var::up));
            row.wTauIDvsJetDown = ratio(n, myCorr->GetTauIDSF_vsJet(kTauID, tau, var::down));
        }
        const float nE = myCorr->GetTauIDSF_vsE(kTauID, tau);
        row.wTauIDvsEUp = ratio(nE, myCorr->GetTauIDSF_vsE(kTauID, tau, var::up));
        row.wTauIDvsEDown = ratio(nE, myCorr->GetTauIDSF_vsE(kTauID, tau, var::down));
        const float nMu = myCorr->GetTauIDSF_vsMu(kTauID, tau);
        row.wTauIDvsMuUp = ratio(nMu, myCorr->GetTauIDSF_vsMu(kTauID, tau, var::up));
        row.wTauIDvsMuDown = ratio(nMu, myCorr->GetTauIDSF_vsMu(kTauID, tau, var::down));
        const float nID = myCorr->GetMuonIDSF(kMuonIDKey, muon);
        row.wMuIDUp = ratio(nID, myCorr->GetMuonIDSF(kMuonIDKey, muon, var::up));
        row.wMuIDDown = ratio(nID, myCorr->GetMuonIDSF(kMuonIDKey, muon, var::down));
        const float nIso = myCorr->GetMuonIDSF(kMuonIsoKey, muon);
        row.wMuIsoUp = ratio(nIso, myCorr->GetMuonIDSF(kMuonIsoKey, muon, var::up));
        row.wMuIsoDown = ratio(nIso, myCorr->GetMuonIDSF(kMuonIsoKey, muon, var::down));
        const auto trig = [&](var v) {
            const float d = myCorr->GetMuonTriggerEff(kMuonTriggerKey, muon.Eta(), muon.Pt(), true, v);
            const float m = myCorr->GetMuonTriggerEff(kMuonTriggerKey, muon.Eta(), muon.Pt(), false, v);
            return m > 0.f ? d / m : 1.f;
        };
        const float nTrig = trig(var::nom);
        row.wMuTrigUp = ratio(nTrig, trig(var::up));
        row.wMuTrigDown = ratio(nTrig, trig(var::down));
        const float nPU = myCorr->GetPUWeight(ev.nTrueInt());
        row.wPUUp = ratio(nPU, myCorr->GetPUWeight(ev.nTrueInt(), var::up));
        row.wPUDown = ratio(nPU, myCorr->GetPUWeight(ev.nTrueInt(), var::down));
        row.wBtagHeavyUp = ratio(btagWeight, bHUp);
        row.wBtagHeavyDown = ratio(btagWeight, bHDown);
        row.wBtagLightUp = ratio(btagWeight, bLUp);
        row.wBtagLightDown = ratio(btagWeight, bLDown);
        row.wISRUp = model_.isrUp;
        row.wISRDown = model_.isrDown;
        row.wFSRUp = model_.fsrUp;
        row.wFSRDown = model_.fsrDown;
        row.wMuRUp = model_.muRUp;
        row.wMuRDown = model_.muRDown;
        row.wMuFUp = model_.muFUp;
        row.wMuFDown = model_.muFDown;
        row.wMuRFUp = model_.muRFUp;
        row.wMuRFDown = model_.muRFDown;
        row.wTopPt = model_.topPt;
        row.wHdampUp = model_.hdampUp;
        row.wHdampDown = model_.hdampDown;
        row.wPdf = model_.pdf;
    }
    row.tauPassMedium = tauMedium;
    row.oppositeSign = os;
    row.nJet = nJet;
    row.nBJet = nBJet;
    row.muonPt = muon.Pt();
    row.muonEta = muon.Eta();
    row.muonPhi = muon.Phi();
    row.tauPt = tau.Pt();
    row.tauMiniAODPt = tau.MiniAODPt();
    row.tauEta = tau.Eta();
    row.tauPhi = tau.Phi();
    row.tauMass = tau.P4().M();
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
    row.metPhi = met.Phi();
    row.jet1Pt = jet1Pt;
    row.lheHT = (!IsDATA && LHE_HT.valid()) ? static_cast<float>(LHE_HT) : -1.f;
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
        const std::vector<GenTau> genTaus = BuildGenTaus(rec, ztt, nullptr);
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
            row.genValidE = r.validE;
            row.genOmegaE = r.omegaE;
            row.genUpsilon = r.upsilon;
            row.genX = r.x;
            row.genVisPt = best->visible.Pt();
            row.genTauPt = best->tau.Pt();
            row.genMatchDR = bestDR;
        }
    }
    (ztt ? recoTausZ_ : (cr ? recoTausCR_ : recoTaus_)).Fill();

    // histogramss  : {OS,SS} x {DM} x {flavour}, nominal rows only.
    if (region < 0)
        return;
    static const char *kCRDir[] = {"", "AntiIso/", "ZeroB/", "AntiIsoZeroB/", "", "AntiIso/"};
    const std::string dir = std::string(cr || ztt ? kCRDir[region] : "") +
                            std::string(os ? "OS" : "SS") + "/DM" + std::to_string(row.tauDM) + "/" +
                            FlavName(row.genFlav) + "/";
    hists.Fill(dir + "Upsilon", row.upsilon, weight, 40, -1.f, 1.f);
    hists.Fill(dir + "TauPt", row.tauPt, weight, 40, 0.f, 200.f);
    hists.Fill(dir + "MVis", row.mVis, weight, 40, 0.f, 200.f);
}

// ZMuMu DY callibration
void TauPolarization::BookZMuMu() {
    zMuMu_ = Output().Book("ZMuMuEvents"); // the histogram directory is ZMuMu
    zMuMu_.Field("weight", zRow_.weight)
        .Field("btagWeight", zRow_.btagWeight)
        .Field("nJet", zRow_.nJet)
        .Field("nBJet", zRow_.nBJet)
        .Field("mll", zRow_.mll)
        .Field("ptll", zRow_.ptll)
        .Field("mu1Pt", zRow_.mu1Pt)
        .Field("mu2Pt", zRow_.mu2Pt)
        .Field("jet1Pt", zRow_.jet1Pt)
        .Field("jet2Pt", zRow_.jet2Pt)
        .Field("ht", zRow_.ht)
        .Field("met", zRow_.met)
        .Field("lheHT", zRow_.lheHT);
}

// two tight isolated opposite-sign muons in the Z window with the mu tau_h jet definition, to check how each DY description models >= 2 reconstructed jets
void TauPolarization::RunZMuMu() {
    auto hists = Hists("ZMuMu");
    const Event ev = GetEvent();
    const float preWeight = IsDATA ? 1.f : MCweight() * ev.GetTriggerLumi(kMuonTrigger) * myCorr->GetPUWeight(ev.nTrueInt());
    if (IsDATA && !myCorr->IsGoldenLumi(RunNumber, ev.lumi()))
        return;
    if (!ev.PassTrigger(kMuonTrigger))
        return;
    const JetViewCollection jets = GetAllJetViews();
    if (!PassNoiseFilterCommon(jets, ev, Event::MET_Type::PUPPI) || !PassJetVetoMap(jets))
        return;

    const MuonViewCollection muons = GetAllMuonViews();
    std::vector<std::size_t> tight;
    int nLooseMuons = 0;
    for (std::size_t i = 0; i < muons.size(); ++i) {
        const auto mu = muons[i];
        if (std::fabs(mu.Eta()) > kMuonEtaMax)
            continue;
        if (mu.Pt() > kVetoLeptonPtMin && mu.PassID(MuonView::MuonID::POG_LOOSE))
            ++nLooseMuons;
        if (mu.Pt() > kSubMuonPtMin && mu.PassID(MuonView::MuonID::POG_TIGHT) && mu.PfRelIso04() < kMuonIsoMax)
            tight.push_back(i);
    }
    if (tight.size() != 2 || nLooseMuons != 2)
        return;
    std::sort(tight.begin(), tight.end(), [&](std::size_t a, std::size_t b) { return muons[a].Pt() > muons[b].Pt(); });
    const auto mu1 = muons[tight[0]], mu2 = muons[tight[1]];
    if (mu1.Pt() < kMuonPtMin || mu1.Charge() * mu2.Charge() >= 0)
        return;
    const ElectronViewCollection electrons = GetAllElectronViews();
    for (std::size_t i = 0; i < electrons.size(); ++i)
        if (electrons[i].Pt() > kVetoLeptonPtMin && electrons[i].PassID(ElectronView::ElectronID::POG_LOOSE))
            return;
    const TrigObjViewCollection trigobjs = GetAllTrigObjViews();
    const auto matched = MatchViewsToTrigObjs(muons, trigobjs, 0.1f, 13);
    const bool m1 = tight[0] < matched.size() && matched[tight[0]] >= 0;
    const bool m2 = tight[1] < matched.size() && matched[tight[1]] >= 0 && mu2.Pt() > kMuonPtMin;
    if (!m1 && !m2)
        return;
    const TLorentzVector ll = mu1.P4() + mu2.P4();
    if (ll.M() < kZMassMin || ll.M() > kZMassMax)
        return;

    const float bWP = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
    int nJet = 0, nBJet = 0;
    float btagWeight = 1.f, ht = 0.f, j1 = 0.f, j2 = 0.f;
    for (std::size_t i = 0; i < jets.size(); ++i) {
        const auto j = jets[i];
        if (j.Pt() < kJetPtMin || std::fabs(j.Eta()) >= kJetEtaMax) 
            continue;
        if (!myCorr->PassJetID(j, JetView::JetID::TIGHT))
            continue;
        if (j.P4().DeltaR(mu1.P4()) < kJetCleanDR || j.P4().DeltaR(mu2.P4()) < kJetCleanDR)
            continue;
        ++nJet;
        ht += j.Pt();
        if (j.Pt() > j1) {
            j2 = j1;
            j1 = j.Pt();
        } else if (j.Pt() > j2) {
            j2 = j.Pt();
        }
        const bool tagged = j.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B) > bWP;
        if (tagged)
            ++nBJet;
        if (!IsDATA) {
            const int flav = static_cast<int>(j.HadronFlavour());
            const float sf = myCorr->GetBTaggingSF(j.Eta(), j.Pt(), flav, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
            if (tagged) {
                btagWeight *= sf;
            } else {
                const float eff = myCorr->GetBTaggingEff(j.Eta(), j.Pt(), flav, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
                if (eff > 0.f && eff < 1.f)
                    btagWeight *= (1.f - sf * eff) / (1.f - eff);
            }
        }
    }

    float weight = preWeight;
    if (!IsDATA) {
        for (const auto &mu : {mu1, mu2}) {
            weight *= myCorr->GetMuonIDSF(kMuonIDKey, mu);
            weight *= myCorr->GetMuonIDSF(kMuonIsoKey, mu);
        }
        // trigger: at least one of the two muons fires IsoMu24
        float effData = 1.f, effMC = 1.f;
        for (const auto &mu : {mu1, mu2}) {
            if (mu.Pt() < kMuonPtMin)
                continue;
            effData *= 1.f - myCorr->GetMuonTriggerEff(kMuonTriggerKey, mu.Eta(), mu.Pt(), true);
            effMC *= 1.f - myCorr->GetMuonTriggerEff(kMuonTriggerKey, mu.Eta(), mu.Pt(), false);
        }
        if (effMC < 1.f)
            weight *= (1.f - effData) / (1.f - effMC);
        weight *= btagWeight;
    }

    ZMuMuRow &row = zRow_;
    row = ZMuMuRow{};
    row.weight = weight;
    row.btagWeight = btagWeight;
    row.nJet = nJet;
    row.nBJet = nBJet;
    row.mll = ll.M();
    row.ptll = ll.Pt();
    row.mu1Pt = mu1.Pt();
    row.mu2Pt = mu2.Pt();
    row.jet1Pt = j1;
    row.jet2Pt = j2;
    row.ht = ht;
    row.met = ev.GetMETVector(Event::MET_Type::PUPPI).Pt();
    row.lheHT = (!IsDATA && LHE_HT.valid()) ? static_cast<float>(LHE_HT) : -1.f;
    zMuMu_.Fill();
    hists.Fill("NJet", static_cast<float>(std::min(nJet, 6)), weight, 7, 0.f, 7.f);
}
