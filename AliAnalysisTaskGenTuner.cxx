/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <Riostream.h>

// ROOT includes
#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TArrayD.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>

// STEER includes
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODDimuon.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskGenTuner.h"


ClassImp(AliAnalysisTaskGenTuner)

//________________________________________________________________________
AliAnalysisTaskGenTuner::AliAnalysisTaskGenTuner() :
AliAnalysisTaskSE(),
fList(0x0),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fMuonTrackCuts(0x0),
fPtCut(-1.),
fWeight(kFALSE),
fDataFile(""),
fPtFunc(0x0),
fPtFuncMC(0x0),
fPtFix(0x0),
fPtFuncNew(0x0),
fPtFixNew(0x0),
fYFunc(0x0),
fYFuncMC(0x0),
fYFix(0x0),
fYFuncNew(0x0),
fYFixNew(0x0),
fcRes(0x0),
fcRat(0x0),
fPtCopyFunc(0x0),
fPtCopyFuncNew(0x0),
fYCopyFunc(0x0),
fYCopyFuncNew(0x0)
{
  /// Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskGenTuner::AliAnalysisTaskGenTuner(const char *name) :
AliAnalysisTaskSE(name),
fList(0x0),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fMuonTrackCuts(0x0),
fPtCut(-1.),
fWeight(kFALSE),
fDataFile(""),
fPtFunc(0x0),
fPtFuncMC(0x0),
fPtFix(0x0),
fPtFuncNew(0x0),
fPtFixNew(0x0),
fYFunc(0x0),
fYFuncMC(0x0),
fYFix(0x0),
fYFuncNew(0x0),
fYFixNew(0x0),
fcRes(0x0),
fcRat(0x0),
fPtCopyFunc(0x0),
fPtCopyFuncNew(0x0),
fYCopyFunc(0x0),
fYCopyFuncNew(0x0)
{
  /// Constructor
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskGenTuner::~AliAnalysisTaskGenTuner()
{
  /// Destructor
  
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fList;
  delete fMuonTrackCuts;
  delete fPtFunc;
  delete fPtFuncMC;
  delete fPtFix;
  delete fPtFuncNew;
  delete fPtFixNew;
  delete fYFunc;
  delete fYFuncMC;
  delete fYFix;
  delete fYFuncNew;
  delete fYFixNew;
  delete fcRes;
  delete fcRat;
  delete fPtCopyFunc;
  delete fPtCopyFuncNew;
  delete fYCopyFunc;
  delete fYCopyFuncNew;
  
}

//___________________________________________________________________________
void AliAnalysisTaskGenTuner::UserCreateOutputObjects()
{
  /// Create histograms
  
  // initialize histos
  fList = new TObjArray(2000);
  fList->SetOwner();
  
  TH1* hPtGen = new TH1D("hPtGen","generated p_{T} distribution;p_{T} (GeV/c);dN/dp_{T}", 300, 0., 30.);
  hPtGen->Sumw2();
  fList->AddAtAndExpand(hPtGen, kPtGen);
  TH1* hPtRec = new TH1D("hPtRec","reconstructed p_{T} distribution;p_{T} (GeV/c);dN/dp_{T}", 300, 0., 30.);
  hPtRec->Sumw2();
  fList->AddAtAndExpand(hPtRec, kPtRec);
  
  TH1* hYGen = new TH1D("hYGen","generated y distribution;y;dN/dy", 250, -4.5, -2.);
  hYGen->Sumw2();
  fList->AddAtAndExpand(hYGen, kYGen);
  TH1* hYRec = new TH1D("hYRec","reconstructed y distribution;y;dN/dy", 250, -4.5, -2.);
  hYRec->Sumw2();
  fList->AddAtAndExpand(hYRec, kYRec);
  
  TH1* hPhiGen = new TH1D("hPhiGen","generated #phi distribution;#phi (deg);dN/d#phi", 360, 0., 360.);
  hPhiGen->Sumw2();
  fList->AddAtAndExpand(hPhiGen, kPhiGen);
  TH1* hPhiRec = new TH1D("hPhiRec","reconstructed #phi distribution;#phi (deg);dN/d#phi", 360, 0., 360.);
  hPhiRec->Sumw2();
  fList->AddAtAndExpand(hPhiRec, kPhiRec);
  
  // recreate the functions for weighting particles because they are lost when streamed to a file (proof and grid mode)
  if (fWeight && fPtFunc && fPtFuncNew && fYFunc && fYFuncNew) {
    fPtCopyFunc = new TF1("fPtCopyFunc", Pt, fPtFunc->GetXmin(), fPtFunc->GetXmax(), fgkNPtParam);
    fPtCopyFunc->SetParameters(fPtFunc->GetParameters());
    fPtCopyFuncNew = new TF1("fPtCopyFuncNew", Pt, fPtFuncNew->GetXmin(), fPtFuncNew->GetXmax(), fgkNPtParam);
    fPtCopyFuncNew->SetParameters(fPtFuncNew->GetParameters());
    fYCopyFunc = new TF1("fYCopyFunc", Y, fYFunc->GetXmin(), fYFunc->GetXmax(), fgkNYParam);
    fYCopyFunc->SetParameters(fYFunc->GetParameters());
    fYCopyFuncNew = new TF1("fYCopyFuncNew", Y, fYFuncNew->GetXmin(), fYFuncNew->GetXmax(), fgkNYParam);
    fYCopyFuncNew->SetParameters(fYFuncNew->GetParameters());
  } else {
    // or disable the weighting
    fWeight = kFALSE;
  }
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::NotifyRun()
{
  /// Prepare processing of new run: load corresponding OADB objects...
  
  // get the trackCuts for this run
  if (!fMuonTrackCuts) AliFatal("You must specify the requested selections (AliMuonTrackCut obj is missing)");
  fMuonTrackCuts->SetRun(fInputHandler);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::UserExec(Option_t *)
{
  /// process event
  
  // get AOD event
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if ( !aod ) return;
  
  // select the centrality range
  Float_t centrality = aod->GetCentrality()->GetCentralityPercentileUnchecked("V0M");
  if (centrality <= fCentMin || centrality > fCentMax) return;
  
  // fill the MC part if running on MC
  TArrayD weight;
  if (MCEvent()) {
    
    weight.Set(MCEvent()->GetNumberOfTracks());
    for (Int_t i = 0; i < MCEvent()->GetNumberOfTracks(); i++) {
      
      AliAODMCParticle *mctrack = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(i));
      
      if (!mctrack->IsPrimary()) continue;
//      if (mctrack->Pt() > 0.9) continue;
//      if (mctrack->Pt() < 0.9 || mctrack->Pt() > 1.) continue;
//      if (mctrack->Y() < -4.4 || (mctrack->Y() > -4.3 && mctrack->Y() < -2.2) || mctrack->Y() > -2.1) continue;
//      if (mctrack->Y() > -4.2 && mctrack->Y() < -2.3) continue;
//      if (mctrack->Y() < -4. || mctrack->Y() > -2.5) continue;
      
      // compute the weights for all primary particles (other weights are 0)
      Double_t y = mctrack->Y();
      Double_t pT = mctrack->Pt();
      if (fWeight) {
	weight[i] = fPtCopyFuncNew->Eval(pT) / fPtCopyFunc->Eval(pT) * fYCopyFuncNew->Eval(y) / fYCopyFunc->Eval(y);
	if (weight[i] < 0.) {
	  AliError(Form("negative weight: y = %g, pT = %g: w = %g", y, pT, weight[i]));
	  weight[i] = 0.;
	}
      } else weight[i] = 1.;
      
      Double_t w = weight[i];
      ((TH1*)fList->UncheckedAt(kPtGen))->Fill(pT, w);
      ((TH1*)fList->UncheckedAt(kYGen))->Fill(y, w);
      ((TH1*)fList->UncheckedAt(kPhiGen))->Fill(mctrack->Phi()*TMath::RadToDeg(), w);
      
    }
    
  }
  
  // fill the reconstructed part
  for (Int_t i = 0; i < aod->GetNumberOfTracks(); i++){
    
    AliAODTrack *track = static_cast<AliAODTrack *>(aod->GetTrack(i));
    
    Double_t pT = track->Pt();
    Int_t mcLabel = track->GetLabel();
    if (!fMuonTrackCuts->IsSelected(track) || (fPtCut > 0. && pT < fPtCut) ||
	(MCEvent() && mcLabel < 0)) continue;
    
//    if (track->Y() > -3.) continue;
    
//    AliAODMCParticle *mctrack = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(mcLabel));
//    if (mctrack->Pt() > 0.9) continue;
//    if (mctrack->Pt() < 0.9 || mctrack->Pt() > 1.) continue;
//    if (mctrack->Y() < -4.4 || (mctrack->Y() > -4.3 && mctrack->Y() < -2.2) || mctrack->Y() > -2.1) continue;
//    if (mctrack->Y() > -4.2 && mctrack->Y() < -2.3) continue;
//    if (mctrack->Y() < -4. || mctrack->Y() > -2.5) continue;
    
    Double_t w = (weight.GetSize() > 0) ? weight[mcLabel] : 1.;
    ((TH1*)fList->UncheckedAt(kPtRec))->Fill(pT, w);
    ((TH1*)fList->UncheckedAt(kYRec))->Fill(track->Y(), w);
    ((TH1*)fList->UncheckedAt(kPhiRec))->Fill(track->Phi()*TMath::RadToDeg(), w);
    
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::Terminate(Option_t *)
{
  /// post-processing
  
  // get current results
  Int_t hIndex[6] = {kPtGen, kYGen, kPhiGen, kPtRec, kYRec, kPhiRec};
  fList = static_cast<TObjArray*>(GetOutputData(1));
  TH1 *h[6];
  for (Int_t i = 0; i < 6; i++) {
    h[i] = static_cast<TH1*>(fList->UncheckedAt(hIndex[i])->Clone());
    h[i]->SetDirectory(0);
  }
  
  // get the fit ranges
  Double_t fitRangeMC[3][2];
  fitRangeMC[0][0] = GetFitLowEdge(*(h[0]));
  fitRangeMC[0][1] = 999.;
  fitRangeMC[1][0] = GetFitLowEdge(*(h[1]));
  fitRangeMC[1][1] = GetFitUpEdge(*(h[1]));
  fitRangeMC[2][0] = h[2]->GetXaxis()->GetXmin();
  fitRangeMC[2][1] = h[2]->GetXaxis()->GetXmax();
  Double_t fitRange[3][2];
  fitRange[0][0] = (fPtCut > 0.) ? TMath::Max(fitRangeMC[0][0], fPtCut) : fitRangeMC[0][0];
  fitRange[0][1] = 15.;
  fitRange[1][0] = -3.98; // not -4. because to the influence of the eta cut
  fitRange[1][1] = -2.5;
  fitRange[2][0] = h[5]->GetXaxis()->GetXmin();
  fitRange[2][1] = h[5]->GetXaxis()->GetXmax();
  
  // compute acc*eff corrections if it is simulated data
  TH1 *hAccEff[3] = {0x0, 0x0, 0x0};
  for (Int_t i = 0; i < 3 && h[i]->GetEntries() > 0; i++) {
    hAccEff[i] = ComputeAccEff(*(h[i]), *(h[i+3]), Form("%sOverGen",h[i+3]->GetName()), "Acc#{times}Eff");
    //hAccEff[i] = static_cast<TH1*>(h[i+3]->Clone(Form("%sOverGen",h[i+3]->GetName())));
    //hAccEff[i]->SetTitle("Acc#{times}Eff");
    //hAccEff[i]->Divide(h[i]);
  }
  
  // get reference data if provided
  TH1 *hRef[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  if (!fDataFile.IsNull()) {
    TFile* dataFile = TFile::Open(fDataFile.Data(),"READ");
    if (!dataFile || !dataFile->IsOpen()) return;
    TObjArray* data = static_cast<TObjArray*>(dataFile->FindObjectAny("Histograms"));
    if (!data) return;
    for (Int_t i = 3; i < 6; i++){
      hRef[i] = static_cast<TH1*>(data->UncheckedAt(hIndex[i])->Clone());
      hRef[i]->SetDirectory(0);
    }
    dataFile->Close();
  }
  
  // compute corrected data
  for (Int_t i = 0; i < 3 && hRef[i+3] && hAccEff[i]; i++) {
    hRef[i] = static_cast<TH1*>(hRef[i+3]->Clone(Form("%sCorr",hRef[i+3]->GetName())));
    hRef[i]->SetTitle("corrected data");
    hRef[i]->Divide(hAccEff[i]);
  }
  
  // normalize histograms
  Bool_t normalized = kFALSE;
  for (Int_t i = 0; i < 3 && hRef[i]; i++) {
    Double_t integral = h[i]->Integral(h[i]->FindBin(fitRange[i][0]), h[i]->FindBin(fitRange[i][1]), "width");
    Double_t norm = (integral != 0.) ? 1./integral : 1.;
    h[i]->Scale(norm);
    h[i+3]->Scale(norm);
    integral = hRef[i]->Integral(hRef[i]->FindBin(fitRange[i][0]), hRef[i]->FindBin(fitRange[i][1]), "width");
    norm = (integral != 0.) ? 1./integral : 1.;
    hRef[i]->Scale(norm);
    hRef[i+3]->Scale(norm);
    normalized = kTRUE;
  }
  
  // compute data/MC ratios
  TH1 *hRat[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  for (Int_t i = 0; i < 6 && hRef[i]; i++) {
    hRat[i] = static_cast<TH1*>(hRef[i]->Clone(Form("%sDataOverMC",hRef[i]->GetName())));
    hRat[i]->SetTitle("data / MC");
    hRat[i]->Divide(h[i]);
  }
  
  // prepare fitting functions
  if (hAccEff[0]) {
    if (fPtFunc) {
      fPtFunc->SetRange(fitRange[0][0], fitRange[0][1]);
      if (fWeight && fPtFuncNew) fPtFunc->SetParameters(fPtFuncNew->GetParameters());
      NormFunc(fPtFunc, fitRange[0][0], fitRange[0][1]);
    }
    if (fPtFuncMC) {
      fPtFuncMC->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
      if (fWeight && fPtFuncNew) fPtFuncMC->SetParameters(fPtFuncNew->GetParameters());
      NormFunc(fPtFuncMC, fitRange[0][0], fitRange[0][1]);
    }
    if (hRef[0] && fPtFuncNew) {
      fPtFuncNew->SetRange(fitRange[0][0], fitRange[0][1]);
      NormFunc(fPtFuncNew, fitRange[0][0], fitRange[0][1]);
    }
    if (!normalized) {
      Double_t integral = h[0]->Integral(h[0]->FindBin(fitRange[0][0]), h[0]->FindBin(fitRange[0][1]), "width");
      if (fPtFunc) fPtFunc->SetParameter(0, fPtFunc->GetParameter(0)*integral);
      if (fPtFuncMC) fPtFuncMC->SetParameter(0, fPtFuncMC->GetParameter(0)*integral);
      if (fPtFuncNew) {
	integral = hRef[0]->Integral(hRef[0]->FindBin(fitRange[0][0]), hRef[0]->FindBin(fitRange[0][1]), "width");
	fPtFuncNew->SetParameter(0, fPtFuncNew->GetParameter(0)*integral);
      }
    }
  }
  if (hAccEff[1]) {
    if (fYFunc) {
      fYFunc->SetRange(fitRange[1][0], fitRange[1][1]);
      if (fWeight && fYFuncNew) fYFunc->SetParameters(fYFuncNew->GetParameters());
      NormFunc(fYFunc, fitRange[1][0], fitRange[1][1]);
    }
    if (fYFuncMC) {
      fYFuncMC->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
      if (fWeight && fYFuncNew) fYFuncMC->SetParameters(fYFuncNew->GetParameters());
      NormFunc(fYFuncMC, fitRange[1][0], fitRange[1][1]);
    }
    if (hRef[1] && fYFuncNew) {
      fYFuncNew->SetRange(fitRange[1][0], fitRange[1][1]);
      NormFunc(fYFuncNew, fitRange[1][0], fitRange[1][1]);
    }
    if (!normalized) {
      Double_t integral = h[1]->Integral(h[1]->FindBin(fitRange[1][0]), h[1]->FindBin(fitRange[1][1]), "width");
      if (fYFunc) fYFunc->SetParameter(0, fYFunc->GetParameter(0)*integral);
      if (fYFuncMC) fYFuncMC->SetParameter(0, fYFuncMC->GetParameter(0)*integral);
      if (fYFuncNew) {
	integral = hRef[1]->Integral(hRef[1]->FindBin(fitRange[1][0]), hRef[1]->FindBin(fitRange[1][1]), "width");
	fYFuncNew->SetParameter(0, fYFuncNew->GetParameter(0)*integral);
      }
    }
  }
  
  // plot results
  fcRes = new TCanvas("cRes", "results", 900, 600);
  fcRes->Divide(3,2);
  fcRes->cd(1);
  gPad->SetLogy();
  if (hAccEff[0]) {
    if (fPtFuncMC) {
      fPtFuncMC->SetLineColor(3);
      fPtFuncMC->SetLineWidth(3);
      h[0]->Fit(fPtFuncMC, "IWLMR", "e0sames");
    } else h[0]->Draw("e0");
    if (fPtFunc) {
      fPtFunc->SetLineColor(4);
      h[0]->Fit(fPtFunc, "IWLMR+");
    }
  }
  if (hRef[0]) {
    hRef[0]->SetLineColor(2);
    if (fPtFuncNew) {
      fPtFuncNew->SetLineColor(2);
      hRef[0]->Fit(fPtFuncNew, "IWLMR", "e0sames");
    } else hRef[0]->Draw("e0sames");
  }
  fcRes->cd(2);
  if (hAccEff[1]) {
    if (fYFuncMC) {
      fYFuncMC->SetLineColor(3);
      fYFuncMC->SetLineWidth(3);
      h[1]->Fit(fYFuncMC, "IWLMR", "e0sames");
    } else h[1]->Draw("e0");
    if (fYFunc) {
      fYFunc->SetLineColor(4);
      h[1]->Fit(fYFunc, "IWLMR+");
    }
  }
  if (hRef[1]) {
    hRef[1]->SetLineColor(2);
    if (fYFuncNew) {
      fYFuncNew->SetLineColor(2);
      hRef[1]->Fit(fYFuncNew, "IWLMR", "e0sames");
    } else hRef[1]->Draw("e0sames");
  }
  for (Int_t i = 2; i < 6; i++) {
    fcRes->cd(i+1);
    if (i == 3) gPad->SetLogy();
    h[i]->Draw("e0");
    if (hRef[i]) {
      hRef[i]->SetLineColor(2);
      hRef[i]->Draw("e0sames");
    }
  }
  
  // normalize functions to their integral in the range used in MC
  if (hAccEff[0] && fPtFunc) {
    fPtFunc->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
    NormFunc(fPtFunc, fitRangeMC[0][0], fitRangeMC[0][1]);
  }
  if (hAccEff[0] && fPtFuncMC) {
    NormFunc(fPtFuncMC, fitRangeMC[0][0], fitRangeMC[0][1]);
  }
  if (hRef[0] && fPtFuncNew) {
    fPtFuncNew->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
    NormFunc(fPtFuncNew, fitRangeMC[0][0], fitRangeMC[0][1]);
  }
  if (hAccEff[1] && fYFunc) {
    fYFunc->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
    NormFunc(fYFunc, fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (hAccEff[1] && fYFuncMC) {
    NormFunc(fYFuncMC, fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (hRef[1] && fYFuncNew) {
    fYFuncNew->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
    NormFunc(fYFuncNew, fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  
  // prepare data/MC function ratios
  TF1 *ptRat = 0x0;
  if (hRat[0] && fPtFunc && fPtFuncNew) {
    ptRat = new TF1("ptRat", PtRat, fitRangeMC[0][0], hRat[0]->GetXaxis()->GetXmax(), 2*fgkNPtParam);
    Double_t p[2*fgkNPtParam];
    fPtFuncNew->GetParameters(p);
    fPtFunc->GetParameters(&(p[fgkNPtParam]));
    ptRat->SetParameters(p);
  }
  TF1 *yRat = 0x0;
  if (hRat[1] && fYFunc && fYFuncNew) {
    yRat = new TF1("yRat", YRat, fitRangeMC[1][0], fitRangeMC[1][1], 2*fgkNYParam);
    Double_t p[2*fgkNYParam];
    fYFuncNew->GetParameters(p);
    fYFunc->GetParameters(&(p[fgkNYParam]));
    yRat->SetParameters(p);
  }
  
  // plot ratios
  fcRat = new TCanvas("cRat", "ratios", 900, 600);
  fcRat->Divide(3,2);
  for (Int_t i = 0; i < 6 && hRat[i]; i++) {
    fcRat->cd(i+1);
    hRat[i]->Draw("e0");
    if (i == 0 && ptRat) ptRat->Draw("sames");
    else if (i == 1 && yRat) yRat->Draw("sames");
  }
  
  // print fitting ranges
  if (fPtFuncMC && fYFuncMC) {
    printf("\npT fitting range MC = [%g, %g]\n", fitRangeMC[0][0], fitRangeMC[0][1]);
    printf("y fitting range MC = [%g, %g]\n\n", fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (fPtFunc && fYFunc) {
    printf("pT fitting range = [%g, %g]\n", fitRange[0][0], fitRange[0][1]);
    printf("y fitting range = [%g, %g]\n\n", fitRange[1][0], fitRange[1][1]);
  }
  
  // print parameters
  Double_t *param = GetOldPtParamMC();
  if (param) {
    printf("Double_t oldPtParamMC[%d] = {", fgkNPtParam);
    for (Int_t i = 0; i < fgkNPtParam-1; i++) printf("%g, ", param[i]);
    printf("%g};\n", param[fgkNPtParam-1]);
  }
  param = GetOldYParamMC();
  if (param) {
    printf("Double_t oldYParamMC[%d] = {", fgkNYParam);
    for (Int_t i = 0; i < fgkNYParam-1; i++) printf("%g, ", param[i]);
    printf("%g};\n\n", param[fgkNYParam-1]);
  }
  param = GetOldPtParam();
  if (param) {
    printf("Double_t oldPtParam[%d] = {", fgkNPtParam);
    for (Int_t i = 0; i < fgkNPtParam-1; i++) printf("%g, ", param[i]);
    printf("%g};\n", param[fgkNPtParam-1]);
  }
  param = GetOldYParam();
  if (param) {
    printf("Double_t oldYParam[%d] = {", fgkNYParam);
    for (Int_t i = 0; i < fgkNYParam-1; i++) printf("%g, ", param[i]);
    printf("%g};\n\n", param[fgkNYParam-1]);
  }
  param = GetNewPtParam();
  if (param) {
    printf("Double_t newPtParam[%d] = {", fgkNPtParam);
    for (Int_t i = 0; i < fgkNPtParam-1; i++) printf("%g, ", param[i]);
    printf("%g};\n", param[fgkNPtParam-1]);
  }
  param = GetNewYParam();
  if (param) {
    printf("Double_t newYParam[%d] = {", fgkNYParam);
    for (Int_t i = 0; i < fgkNYParam-1; i++) printf("%g, ", param[i]);
    printf("%g};\n\n", param[fgkNYParam-1]);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::SetPtParam(const Double_t *pOld, const Bool_t *fixOld, const Double_t *pNew,
					 const Bool_t *fixNew, Double_t min, Double_t max)
{
  /// create the function(s) to fit the generated pT distribution with the given parameters
  /// if the new parameters are given, the new function is used to weight the tracks
  
  // prepare the old fitting functions
  if (!fPtFunc) fPtFunc = new TF1("fPtFunc", Pt, min, max, fgkNPtParam);
  else fPtFunc->SetRange(min, max);
  if (!fPtFuncMC) fPtFuncMC = new TF1("fPtFuncMC", Pt, min, max, fgkNPtParam);
  else fPtFuncMC->SetRange(min, max);
  
  // set the current parameters
  if (pOld) {
    fPtFunc->SetParameters(pOld);
    fPtFuncMC->SetParameters(pOld);
  }
  
  // fix some of them if required
  if (fixOld) {
    if (!fPtFix) fPtFix = new Bool_t[fgkNPtParam];
    for (Int_t i = 0; i < fgkNPtParam; i++) {
      fPtFix[i] = fixOld[i];
      if (fPtFix[i]) {
	fPtFunc->FixParameter(i, fPtFunc->GetParameter(i));
	fPtFuncMC->FixParameter(i, fPtFuncMC->GetParameter(i));
      }
    }
  }
  
  // normalize the functions
  NormFunc(fPtFunc, min, max);
  NormFunc(fPtFuncMC, min, max);
  
  if (pNew) {
    
    // prepare the new fitting function if needed
    if (!fPtFuncNew) fPtFuncNew = new TF1("fPtFuncNew", Pt, min, max, fgkNPtParam);
    else fPtFuncNew->SetRange(min, max);
    
    // set the current parameters
    fPtFuncNew->SetParameters(pNew);
    
    // fix some of them if required
    if (fixNew) {
      if (!fPtFixNew) fPtFixNew = new Bool_t[fgkNPtParam];
      for (Int_t i = 0; i < fgkNPtParam; i++) {
	fPtFixNew[i] = fixNew[i];
	if (fPtFixNew[i]) fPtFuncNew->FixParameter(i, fPtFuncNew->GetParameter(i));
      }
    }
    
    // normalize the function
    NormFunc(fPtFuncNew, min, max);
    
  } else {
    
    // or make sure the new function does not exist if not required
    delete fPtFuncNew;
    fPtFuncNew = 0x0;
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::SetYParam(const Double_t *pOld, const Bool_t *fixOld, const Double_t *pNew,
					const Bool_t *fixNew, Double_t min, Double_t max)
{
  /// create the function(s) to fit the generated y distribution with the given parameters
  /// if the new parameters are given, the new function is used to weight the tracks
  
  // prepare the old fitting functions
  if (!fYFunc) fYFunc = new TF1("fYFunc", Y, min, max, fgkNYParam);
  else fYFunc->SetRange(min, max);
  if (!fYFuncMC) fYFuncMC = new TF1("fYFuncMC", Y, min, max, fgkNYParam);
  else fYFuncMC->SetRange(min, max);
  
  // set the current parameters
  if (pOld) {
    fYFunc->SetParameters(pOld);
    fYFuncMC->SetParameters(pOld);
  }
  
  // fix some of them if required
  if (fixOld) {
    if (!fYFix) fYFix = new Bool_t[fgkNYParam];
    for (Int_t i = 0; i < fgkNYParam; i++) {
      fYFix[i] = fixOld[i];
      if (fYFix[i]) {
	fYFunc->FixParameter(i, fYFunc->GetParameter(i));
	fYFuncMC->FixParameter(i, fYFuncMC->GetParameter(i));
      }
    }
  }
  
  // normalize the function
  NormFunc(fYFunc, min, max);
  NormFunc(fYFuncMC, min, max);
  
  if (pNew) {
    
    // prepare the new fitting function if needed
    if (!fYFuncNew) fYFuncNew = new TF1("fYFuncNew", Y, min, max, fgkNYParam);
    else fYFuncNew->SetRange(min, max);
    
    // set the current parameters
    fYFuncNew->SetParameters(pNew);
    
    // fix some of them if required
    if (fixNew) {
      if (!fYFixNew) fYFixNew = new Bool_t[fgkNYParam];
      for (Int_t i = 0; i < fgkNYParam; i++) {
	fYFixNew[i] = fixNew[i];
	if (fYFixNew[i]) fYFuncNew->FixParameter(i, fYFuncNew->GetParameter(i));
      }
    }
    
    // normalize the function
    NormFunc(fYFuncNew, min, max);
    
  } else {
    
    // or make sure the new function does not exist if not required
    delete fYFuncNew;
    fYFuncNew = 0x0;
    
  }
  
}

//________________________________________________________________________
TH1* AliAnalysisTaskGenTuner::ComputeAccEff(TH1 &hGen, TH1 &hRec, const Char_t *name, const Char_t *title)
{
  /// Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  /// Result is identical to divide histograms with option "B", except here error is forced > 1/gen
  
  Int_t nbins = hGen.GetNbinsX();
  TH1* hAcc = new TH1D(name,title, nbins, hGen.GetXaxis()->GetXmin(), hGen.GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t accEff = 0.;
    Double_t accEffErr = 0.;
    Double_t gen = hGen.GetBinContent(i);
    if (gen > 0.) {
      Double_t rec = hRec.GetBinContent(i);
      Double_t genErr = hGen.GetBinError(i);
      Double_t recErr = hRec.GetBinError(i);
      accEff = rec/gen;
      //Double_t accEffErr2 = ((1.-2.*accEff)*recErr*recErr + accEff*accEff*genErr*genErr)/(gen*gen);
      //accEffErr = TMath::Max(accEff*genErr/gen, TMath::Sqrt(TMath::Abs(accEffErr2)));
      accEffErr = TMath::Max(recErr/gen, accEff*genErr/gen);
    }
    hAcc->SetBinContent(i, accEff);
    hAcc->SetBinError(i, accEffErr);
  }
  
  return hAcc;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::Pt(const Double_t *x, const Double_t *p)
{
  /// generated pT fit function
  Double_t pT = *x;
  return p[0] * (1. / TMath::Power(p[1] + TMath::Power(pT,p[2]), p[3]) + p[4] * TMath::Exp(p[5]*pT));
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::Y(const Double_t *x, const Double_t *p)
{
  /// generated y fit function
  Double_t y = *x;
  Double_t arg = y/p[7];
  return p[0] * (p[1] * (1. + p[2]*y + p[3]*y*y + p[4]*y*y*y + p[5]*y*y*y*y) + p[6]*TMath::Exp(-0.5*arg*arg));
}

//________________________________________________________________________

Double_t AliAnalysisTaskGenTuner::PtRat(const Double_t *x, const Double_t *p)
{
  /// generated pT fit function ratio
  const Double_t *p1 = &(p[0]);
  const Double_t *p2 = &(p[fgkNPtParam]);
  return Pt(x,p1) / Pt(x,p2);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::YRat(const Double_t *x, const Double_t *p)
{
  /// generated y fit function ratio
  const Double_t *p1 = &(p[0]);
  const Double_t *p2 = &(p[fgkNYParam]);
  return Y(x,p1) / Y(x,p2);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::GetFitLowEdge(TH1 &h)
{
  /// adjust the lower edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindFirstBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0+1)) binAbove0++;
  return h.GetBinLowEdge(binAbove0);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::GetFitUpEdge(TH1 &h)
{
  /// adjust the upper edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindLastBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0-1)) binAbove0--;
  return h.GetBinLowEdge(binAbove0+1);
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::NormFunc(TF1 *f, Double_t min, Double_t max)
{
  /// normalize the function to its integral in the given range
  Double_t integral = f->Integral(min, max);
  if (integral != 0.) f->SetParameter(0, f->GetParameter(0)/integral);
}

