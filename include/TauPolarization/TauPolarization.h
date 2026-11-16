#ifndef TauPolarization_TauPolarization_h
#define TauPolarization_TauPolarization_h

#include <cstddef>
#include <vector>

#include <AnalyzerFramework/AnalyzerCore.h>

#include "MyCorrection.h"

// tau polarization in W -> tau nu. Tasks (--userflags):
//   GenStudy : generator-level polarimeter (MC only)
//   MuTauSR  : reco-level ttbar -> mu tau_h selection (data and MC).
class TauPolarization final : public AnalyzerCore {
public:
    void initializeAnalyzer() override;
    void executeEvent() override;

private:
    void RegisterTasks();

    void ValidateGenStudy();
    void BookGenStudy();
    void RunGenStudy();

    void BookMuTauSR();
    void RunMuTauSR();

    struct GenTauRow {
        float weight = 0.f;
        int charge = 0;
        int origin = 0;          // 0 = W, 1 = t -> W
        int ancestorPdgId = 0;
        int genDM = 0;           // from GenPart daughters (-1 leptonic, 15 other)
        int visTauDM = -1;       // GenVisTau_status of the matched GenVisTau, -1 if none
        int chargedPdgId = 0;
        int nPhotons = 0;
        bool valid = false;      // omega computed
        bool clamped = false;    // |omega| > 1 from stored-precision kinematics, clamped
        float omega = -9.f;
        float cosThetaStar = -9.f;
        float cosPsi = -9.f;
        float x = -9.f;
        float upsilon = -9.f;
        float tauPt = 0.f, tauEta = 0.f, tauE = 0.f;
        float visPt = 0.f, visEta = 0.f, visE = 0.f, visM = 0.f;
        float visTauDR = -1.f;   // dR(daughter sum, matched GenVisTau)
        float visTauRelDE = -9.f;// (E_daughters - E_GenVisTau) / E_GenVisTau
        float nuResidual = -1.f; // |(tau - vis) - stored nu_tau| / E_tau, -1 if nu not stored
    };

    // One row per selected mu tau_h event (OS and SS), before the sign cut.
    struct RecoTauRow {
        float weight = 0.f;      // full central weight (1 on data)
        bool oppositeSign = false;
        int nJet = 0, nBJet = 0;
        float muonPt = 0.f, muonEta = 0.f;
        float tauPt = 0.f;       // TES-corrected
        float tauMiniAODPt = 0.f;
        float tauEta = 0.f;
        int tauCharge = 0;
        int tauDM = -1;
        int tauDMUParT = -1;
        float leadTkRatio = -1.f;// Tau_leadTkPtOverTauPt (relative to the MiniAOD pT)
        float upsilon = -9.f;    // 2 R - 1 with R rescaled to the TES-corrected pT
        float upsilonTesUp = -9.f, upsilonTesDown = -9.f;
        float tauPtTesUp = 0.f, tauPtTesDown = 0.f;
        float mVis = 0.f, mT = 0.f, met = 0.f;
        int genFlav = -1;        // Tau_genPartFlav, -1 on data
        // Matched generator tau from a W (genFlav 5 only), -9/-1 otherwise.
        int genDM = -1;
        int genCharge = 0;
        int genOrigin = -1;
        int genNPhotons = -1;
        bool genValid = false;
        bool genClamped = false;
        float genOmega = -9.f;
        float genUpsilon = -9.f;
        float genX = -9.f;
        float genVisPt = 0.f;
        float genTauPt = 0.f;
        float genMatchDR = -1.f;
    };

    GenTauRow genRow_;
    RNTupleHandle genTaus_;
    RecoTauRow recoRow_;
    RNTupleHandle recoTaus_;
};

#endif
