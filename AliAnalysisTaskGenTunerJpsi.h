#ifndef ALIANALYSISTASKGENTUNERJPSI_H
#define ALIANALYSISTASKGENTUNERJPSI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskGenTunerJpsi
/// \brief task to tune the muon or JPsi pt/y generated distributions
//Author: Philippe Pillot - SUBATECH Nantes

#include <TObject.h>
#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"

class TH1;
class TF1;
class TCanvas;

//________________________________________________________________________
class AliAnalysisTaskGenTunerJpsi : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskGenTunerJpsi();
  AliAnalysisTaskGenTunerJpsi(const char *name);
  virtual ~AliAnalysisTaskGenTunerJpsi();
  
  virtual void   UserCreateOutputObjects();
  virtual void   NotifyRun();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  /// Select tracks in the given centrality range
  void SelectCentrality(Double_t min, Double_t max) {fCentMin = min; fCentMax = max;}
  
  // set standard cuts to select tracks to be considered
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  
  /// set the muon low pT cut
  void SetMuonPtCut(Double_t cut) {fPtCut = cut;}
  
  /// weight simulated/reconstructed particles using using given functions
  void Weight(Bool_t flag) {fWeight = flag;}
  
  // set the pt ref. histo to reweight MC
  void SetPtRefHisto(TH1* hptRef) {if(hptRef) fHptRef = static_cast<TH1*>(hptRef->Clone());}

  // set the Y ref. histo to reweight MC
  void SetYRefHisto(TH1* hyRef) {if(hyRef) fHyRef = static_cast<TH1*>(hyRef->Clone());}

  // set pt binning
  void SetPtBin(Int_t nofbin,Double_t* bin); 

  // set y binning
  void SetYBin(Int_t nofbin,Double_t* bin ); 
  
  // create the original function with the parameters used in simulation to generate the pT distribution
  void SetOriginPtFunc(TString formula, const Double_t *param, const Bool_t *fixParam, Double_t xMin, Double_t xMax);
  // create the new function with its initial parameters to fit the generated/weighted pT distribution
  void SetNewPtFunc(TString formula, const Double_t *param, const Bool_t *fixParam, Double_t xMin, Double_t xMax);
  
  // create the original function with the parameters used in simulation to generate the y distribution
  void SetOriginYFunc(TString formula, const Double_t *param, const Bool_t *fixParam, Double_t xMin, Double_t xMax);
  // create the new function with its initial parameters to fit the generated/weighted y distribution
  void SetNewYFunc(TString formula, const Double_t *param, const Bool_t *fixParam, Double_t xMin, Double_t xMax);
 
  /// get the generated pT fit function with current parameters
  TF1* GetCurrentPtFunc() {return fPtFunc;}
  /// get the generated pT fit function with current parameters fitted in the MC range
  TF1* GetCurrentPtFuncMC() {return fPtFuncMC;}
  /// get the generated pT fit function with new parameters
  TF1* GetNewPtFunc() {return fPtFuncNew;}
  
  /// get the generated y fit function with current parameters
  TF1* GetCurrentYFunc() {return fYFunc;}
  /// get the generated y fit function with current parameters fitted in the MC range
  TF1* GetCurrentYFuncMC() {return fYFuncMC;}
  /// get the generated y fit function with new parameters
  TF1* GetNewYFunc() {return fYFuncNew;}
  
  // get canvas containing generated and reconstructed distributions
  TCanvas* GetResults() {return fcRes;}
  // get canvas containing data/MC ratios
  TCanvas* GetRatios() {return fcRat;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskGenTunerJpsi(const AliAnalysisTaskGenTunerJpsi& rhs);
  /// Not implemented
  AliAnalysisTaskGenTunerJpsi& operator = (const AliAnalysisTaskGenTunerJpsi& rhs);
  
  // Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  TH1* ComputeAccEff(TH1 &hGen, TH1 &hRec, const Char_t *name, const Char_t *title);
  
  // adjust the lower edge of the fit range according to the content of the histogram
  Double_t GetFitLowEdge(TH1 &h);
  // adjust the upper edge of the fit range according to the content of the histogram
  Double_t GetFitUpEdge(TH1 &h);
  
  // normalize the function to its integral in the given range
  void NormFunc(TF1 *f, Double_t min, Double_t max);
  
  // generated pT fit function ratio
  Double_t PtRat(const Double_t *x, const Double_t *p);
  // generated y fit function ratio
  Double_t YRat(const Double_t *x, const Double_t *p);
  
private:
  
   enum eList {
    kPtGen   = 0, ///< pT distribution of generated particle
    kPtRec   = 1, ///< pT distribution of reconstructed particle
    kYGen    = 2, ///< y distribution of generated particle
    kYRec    = 3,  ///< y distribution of reconstructed particle
    // kPhiGen  = 4, ///< phi distribution of generated particle
    // kPhiRec  = 5  ///< phi distribution of reconstructed particle
  };
  
  TObjArray*  fList; //!< List of output object
  
  Double_t  fCentMin;               ///< select centrality > fCentMin
  Double_t  fCentMax;               ///< select centrality <= fCentMax
  
  AliMuonTrackCuts* fMuonTrackCuts; ///< cuts to select tracks to be considered
  Double_t fPtCut;                  ///< muon low pT cut
  Double_t fGenPtCut;               ///< generated muon low pT cut
  Bool_t   fWeight;                 ///< weight simulated/reconstructed particles using using given functions
  // TString  fDataFile;               ///< data file used in terminate to tune the generated distributions
  TF1     *fPtFuncOld;              ///< original generated pT function with original parameters
  TF1     *fPtFuncNew;              ///< new generated pT fit function with new parameters
  TF1     *fPtFunc;                 //!< current generated pT fit function with current parameters
  TF1     *fPtFuncMC;               //!< current generated pT fit function with current parameters in MC range
  TF1     *fYFuncOld;               ///< original generated y function with original parameters
  TF1     *fYFuncNew;               ///< new generated y fit function with new parameters
  TF1     *fYFunc;                  //!< current generated y fit function with current parameters
  TF1     *fYFuncMC;                //!< current generated y fit function with current parameters in MC range
  TCanvas *fcRes;                   //!< generated and reconstructed distributions
  TCanvas *fcRat;                   //!< data/MC ratios

  //______ Specific datamembers for GenTunerJpsi compare to Gentuner _____
  Int_t     fPtNofBin;              // Number of Pt bin
  Int_t     fYNofBin;               // Number of Y bin 
  TH1*      fHptRef;                ///< pt ref. data histo        
  TH1*      fHyRef;                 ///< y ref. data histo
  Double_t* fPtBin;                 //[fPtNofBin]
  Double_t* fYBin;                  //[fYNofBin] 
 
  
  
  ClassDef(AliAnalysisTaskGenTunerJpsi, 1);
};

//________________________________________________________________________
inline void AliAnalysisTaskGenTunerJpsi::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set standard cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}

#endif

