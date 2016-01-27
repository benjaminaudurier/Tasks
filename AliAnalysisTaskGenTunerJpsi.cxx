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
#include <TLegend.h>

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
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuSpectra.h"

#include "AliAnalysisTaskGenTunerJpsi.h"

ClassImp(AliAnalysisTaskGenTunerJpsi)

//________________________________________________________________________
AliAnalysisTaskGenTunerJpsi::AliAnalysisTaskGenTunerJpsi() :
AliAnalysisTaskSE(),
fList(0x0),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fMuonTrackCuts(0x0),
fPtCut(-1.),
fGenPtCut(-1.),
fWeight(kFALSE),
fPtFuncOld(0x0),
fPtFuncNew(0x0),
fPtFunc(0x0),
fPtFuncMC(0x0),
fYFuncOld(0x0),
fYFuncNew(0x0),
fYFunc(0x0),
fYFuncMC(0x0),
fcRes(0x0),
fcRat(0x0),
fPtNofBin(0),
fYNofBin(0),
fPtBin(0),
fYBin(0),
fHptRef(0x0),
fHyRef(0x0)
{
  /// Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskGenTunerJpsi::AliAnalysisTaskGenTunerJpsi(const char *name) :
AliAnalysisTaskSE(name),
fList(0x0),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fMuonTrackCuts(0x0),
fPtCut(-1.),
fGenPtCut(-1.),
fWeight(kFALSE),
fPtFuncOld(0x0),
fPtFuncNew(0x0),
fPtFunc(0x0),
fPtFuncMC(0x0),
fYFuncOld(0x0),
fYFuncNew(0x0),
fYFunc(0x0),
fYFuncMC(0x0),
fcRes(0x0),
fcRat(0x0),
fPtNofBin(0),
fPtBin(0),
fYNofBin(0),
fYBin(0),                  
fHptRef(0x0),
fHyRef(0x0)
{
  /// Constructor
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskGenTunerJpsi::~AliAnalysisTaskGenTunerJpsi()
{
  /// Destructor
  
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fList;
  delete fMuonTrackCuts;
  delete fPtFuncOld;
  delete fPtFuncNew;
  delete fPtFunc;
  delete fPtFuncMC;
  delete fYFuncOld;
  delete fYFuncNew;
  delete fYFunc;
  delete fYFuncMC;
  delete fcRes;
  delete fcRat;
  delete fHptRef;
  delete fHyRef;
  delete[] fPtBin;                 
  delete[] fYBin;                 
}

//___________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::UserCreateOutputObjects()
{
  /// Create histograms
  
  // initialize histos
  fList = new TObjArray(2000);
  fList->SetOwner();
  
Int_t nofptbin = fPtNofBin-1;
Int_t nofybin = fYNofBin-1;

if (nofptbin==0 || nofybin==0 )
{
  AliError("I don't have the number of bin ...");
  return;
}

  TH1* hPtGen = new TH1D("hPtGen","generated p_{T} distribution;p_{T} (GeV/c);dN/dp_{T}",nofptbin,fPtBin);
  hPtGen->Sumw2();// Create structure to store sum of squares of weights
  fList->AddAtAndExpand(hPtGen, kPtGen);// Double array size if reach the end
  TH1* hPtRec = new TH1D("hPtRec","reconstructed p_{T} distribution;p_{T} (GeV/c);dN/dp_{T}",nofptbin,fPtBin);
  hPtRec->Sumw2();// Create structure to store sum of squares of weights
  fList->AddAtAndExpand(hPtRec, kPtRec);// Double array size if reach the end
  
  TH1* hYGen = new TH1D("hYGen","generated y distribution;y;dN/dy",nofybin,fYBin);
  hYGen->Sumw2();
  fList->AddAtAndExpand(hYGen, kYGen);
  TH1* hYRec = new TH1D("hYRec","reconstructed y distribution;y;dN/dy",nofybin,fYBin);
  hYRec->Sumw2();
  fList->AddAtAndExpand(hYRec, kYRec);

  if(fHptRef)fHptRef->Sumw2();// Set histo error
  if(fHyRef)fHyRef->Sumw2();// set histo error

  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::NotifyRun()
{
  /// Prepare processing of new run: load corresponding OADB objects...
  
  // get the trackCuts for this run
  if (!fMuonTrackCuts) AliFatal("You must specify the requested selections (AliMuonTrackCut obj is missing)");
  fMuonTrackCuts->SetRun(fInputHandler);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::UserExec(Option_t *)
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
  
  if (!MCEvent()) return; 

    
  weight.Set(MCEvent()->GetNumberOfTracks());

  //__________Generated tracks part
  for (Int_t i = 0; i < MCEvent()->GetNumberOfTracks(); i++) {
    
    AliAODMCParticle *mctrack = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(i));
    Int_t mctrackPID = mctrack->PdgCode();
    
    if (!mctrack->IsPrimary()) continue; // Select primary
    if(mctrackPID!=443) continue; // Select J/psi

    
    // compute the weights for all primary particles (other weights are 0)
    Double_t y = mctrack->Y();
    Double_t pT = mctrack->Pt();

    // Cut on j/psi
    if( pT > 8. || pT < 0. || y < -4. || y > -2.5 ) continue; 
    
    if (fWeight && fPtFuncOld && fPtFuncNew && fYFuncOld && fYFuncNew) 
    {
      weight[i] = fPtFuncNew->Eval(pT) / fPtFuncOld->Eval(pT) * fYFuncNew->Eval(y) / fYFuncOld->Eval(y);       
      if (weight[i] < 0.) 
      {
        AliError(Form("negative weight: y = %g, pT = %g: w = %g", y, pT, weight[i]));
        weight[i] = 0.;
      }
    } 
    else weight[i] = 1.;
    
    Double_t w = weight[i];

    ((TH1*)fList->UncheckedAt(kPtGen))->Fill(pT, w);
    ((TH1*)fList->UncheckedAt(kYGen))->Fill(y, w); // fill YGenHisto
  }
    
   

  //__________Rec. Part
  for (Int_t i = 0; i < aod->GetNumberOfTracks(); i++)
  {
    AliAODTrack *track1 = static_cast<AliAODTrack*>(aod->GetTrack(i)); // track 1
    
    Int_t mcLabel1 = track1->GetLabel();
    if (!fMuonTrackCuts->IsSelected(track1) || (MCEvent() && mcLabel1 < 0)) continue;
    //    Physic selection                      track associated with mc

    AliAODMCParticle *mctrack = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(mcLabel1)); // MC info associated with track1
    if((mctrack->PdgCode()!=13 && mctrack->PdgCode()!=-13)) continue; // Cut on muons

    //Check mother
    Int_t mother = mctrack->GetMother(); // mother label
    if(mother<0) continue;
    AliAODMCParticle *mcmothertrack = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(mother)); // MC mother trac
    if(mcmothertrack->PdgCode()!=443 ) continue; // Check if mother is a J/psi

      // Get 4-Vector  
      TLorentzVector vec(track1->Px(),track1->Py(),track1->Pz(),track1->E());
      // cout << "vec.M() = " <<vec.M()<< endl;
    
    for (Int_t j = i+1; j < aod->GetNumberOfTracks(); j++)// Loop over second paricles 
    {
      AliAODTrack *track2 = static_cast<AliAODTrack*>(aod->GetTrack(j)); // track 2
    
      Int_t mcLabel2 = track2->GetLabel();
      if (!fMuonTrackCuts->IsSelected(track2) || (MCEvent() && mcLabel2 < 0) || (track2->Charge()*track1->Charge()>=0)   ) continue;
      //    Physic selection                      track associated with mc part.    opposite sign  
      
      // MC info associated with track2    
      AliAODMCParticle *mctrack2 = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(mcLabel2)); 
      
      if((mctrack2->PdgCode()!=13 && mctrack2->PdgCode()!=-13)  ) continue; 
      //                  Cut muons                               

      //Check mother
      Int_t mother2 = mctrack2->GetMother(); // mother label
      if(mother2<0) continue;
      AliAODMCParticle *mcmothertrack2 = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(mother2)); // MC mother trac
      if(mcmothertrack2->PdgCode()!=443 ) continue; // Check if mother is a J/psi

      // Get 4-Vector
      TLorentzVector vec2(track2->Px(),track2->Py(),track2->Pz(),track2->E());

      // Add the two 4-vectors.
      TLorentzVector Dimu = vec2+vec;

      // Cut on paire
      if( Dimu.Pt() > 8. || Dimu.Pt()<0. || Dimu.Rapidity() < -4. || Dimu.Rapidity()>-2.5 ) continue; 

      if (fWeight && fPtFuncOld && fPtFuncNew && fYFuncOld && fYFuncNew) 
      {
        weight[i] = fPtFuncNew->Eval(Dimu.Pt()) / fPtFuncOld->Eval(Dimu.Pt()) * fYFuncNew->Eval(Dimu.Rapidity()) / fYFuncOld->Eval(Dimu.Rapidity());        
        if (weight[i] < 0.) 
        {
          AliError(Form("negative weight: Dimu.Rapidity() = %g, Dimu.Pt() = %g: w = %g", Dimu.Rapidity(), Dimu.Pt(), weight[i]));
          weight[i] = 0.;
        }
      } 
      else weight[i] = 1.;
      //________
      Double_t w = weight[i];

      ((TH1*)fList->UncheckedAt(kPtRec))->Fill(Dimu.Pt(), w);
      ((TH1*)fList->UncheckedAt(kYRec))->Fill(Dimu.Rapidity(), w); // fill YGenHisto
    } 
  } 
  //__________

   
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::Terminate(Option_t *)
{
  /// post-processing
  
  //__________get current results
  Int_t hIndex[4] = {kPtGen, kYGen, /*kPhiGen,*/ kPtRec, kYRec/*, kPhiRec*/};
  fList = static_cast<TObjArray*>(GetOutputData(1));
  TH1 *h[4]; // List to store histo
  for (Int_t i = 0; i < 4; i++) 
  {
    h[i] = static_cast<TH1*>(fList->UncheckedAt(hIndex[i])->Clone());
    h[i]->SetDirectory(0);
    h[i]->Scale(1,"width");// Normalize
  }

  //__________get the fit ranges
  Double_t fitRangeMC[2][2];
  fitRangeMC[0][0] = GetFitLowEdge(*(h[0]));
  fitRangeMC[0][1] = 50.;
  fitRangeMC[1][0] = GetFitLowEdge(*(h[1]));
  fitRangeMC[1][1] = GetFitUpEdge(*(h[1]));
  Double_t fitRange[2][2];
  fitRange[0][0] = (fPtCut > 0.) ? TMath::Max(fitRangeMC[0][0], fPtCut) : fitRangeMC[0][0];
  fitRange[0][1] = 8.;
  fitRange[1][0] = -4.; // not -4. because to the influence of the eta cut
  fitRange[1][1] = -2.5;
  //__________
  //
  //__________compute acc*eff corrections if it is simulated data
  TH1 *hAccEff[2] = {0x0,       0x0/*,      0x0*/};
  //                 AccEffPt   AccEffY   AccEffPhi

  for (Int_t i = 0; i < 2 && h[i]->GetEntries() > 0; i++) 
  {
    hAccEff[i] = ComputeAccEff(*(h[i]), *(h[i+2]), Form("%sOverGen",h[i+2]->GetName()), "Acc#{times}Eff");
    hAccEff[i]->SetTitle("Acc#{times}Eff");
    fList->Add(hAccEff[i]);
  }
  //__________
  
  //__________get reference data if provided
  TH1 *hRef[4] = {0x0,            0x0,             0x0,          0x0};
  //              ptAccEffCorr.   yAccEffCorr.     RefPtHisto,   RefYHisto, 

  // J/psi ref. histo
  if(fHptRef) hRef[2] = static_cast<TH1*>(fHptRef->Clone());
  if(fHyRef) hRef[3]  = static_cast<TH1*>(fHyRef->Clone());
  //__________

  //__________compute corrected data
  for (Int_t i = 0; i < 2 && hRef[i+2] && hAccEff[i]; i++) 
  {
    hRef[i] = static_cast<TH1*>(hRef[i+2]->Clone());
    
    if(!hRef[i])
    {
      cout << "Can't clone ref. histo " << i << endl;
      continue;
    }
    hRef[i]->SetTitle("corrected data");
    if(!hRef[i]->Divide(hAccEff[i]))
    {
      cout << "Can't divide " << i << endl;
      continue;
    }
  }
  //__________
  
  //__________normalize histograms
  Bool_t normalized = kFALSE;
  for (Int_t i = 0; i < 2 && hRef[i]; i++) 
  {
    Double_t integral = h[i]->Integral("width");
    Double_t norm = (integral != 0.) ? 1./integral : 1.;
    h[i]->Scale(norm);
    h[i+2]->Scale(norm);
    integral = hRef[i]->Integral("width");
    norm = (integral != 0.) ? 1./integral : 1.;
    hRef[i]->Scale(norm);
    hRef[i+2]->Scale(norm);
    normalized = kTRUE; 
  }
  //_________

  //__________compute dataCorr/MC ratios
  TH1 *hRat[4] = {0x0,                  0x0,                 0x0,              0x0            };
  //              ptAccEffCorr./ptGen   yAccEffCorr./yGen    RefPtHisto/RecPt  RefYHisto/RecY 
  for (Int_t i = 0; i < 4 && hRef[i]; i++) 
  {
    hRat[i] = static_cast<TH1*>(hRef[i]->Clone());
    hRat[i]->SetTitle("data / MC");
    hRat[i]->Divide(h[i]);
    if(!hRat[i]) cout << Form("Cannot divide histo %d",i) << endl;
    ;
  }
  //__________

  //__________prepare fitting functions
  if (hAccEff[0]) // Pt AccEff histo
  {
    cout << "prepare fitting functions for Pt AccEff histo... " << endl;

    if (fWeight && fPtFuncNew) {
      fPtFunc = new TF1(*fPtFuncNew);
      fPtFuncMC = new TF1(*fPtFuncNew);
    } else if (!fWeight && fPtFuncOld) {
      fPtFunc = new TF1(*fPtFuncOld);
      fPtFuncMC = new TF1(*fPtFuncOld);
    }
    if (fPtFunc) {
      fPtFunc->SetName("fPtFunc");
      fPtFunc->SetRange(fitRange[0][0], fitRange[0][1]);
      NormFunc(fPtFunc, fitRange[0][0], fitRange[0][1]);
    }
    if (fPtFuncMC) {
      fPtFuncMC->SetName("fPtFuncMC");
      fPtFuncMC->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
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
  
  if (hAccEff[1]) // Y AccEff histo
  {
    cout << "prepare fitting functions for Y AccEff histo..." << endl;

   if (fWeight && fYFuncNew) {
      fYFunc = new TF1(*fYFuncNew);
      fYFuncMC = new TF1(*fYFuncNew);
    } else if (!fWeight && fYFuncOld) {
      fYFunc = new TF1(*fYFuncOld);
      fYFuncMC = new TF1(*fYFuncOld);
    }
    if (fYFunc) {
      fYFunc->SetName("fYFunc");
      fYFunc->SetRange(fitRange[1][0], fitRange[1][1]);
      NormFunc(fYFunc, fitRange[1][0], fitRange[1][1]);
    }
    if (fYFuncMC) {
      fYFuncMC->SetName("fYFuncMC");
      fYFuncMC->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
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
  //__________
  
  
  //__________plot results
  fcRes = new TCanvas("cRes", "results", 900, 600);
  fcRes->Divide(2,2);

  //------Pad 1
  fcRes->cd(1);
  gPad->SetLogy();
  TLegend*leg = new TLegend(0.48,0.7,0.70,0.9);

  if (hAccEff[0]) 
  {
    if (fPtFuncMC) 
    {
      fPtFuncMC->SetLineColor(3);// Green
      fPtFuncMC->SetLineWidth(3);
      h[0]->Fit(fPtFuncMC, "IWLR", "e0sames");// Gen. pt histo.
      leg->AddEntry(fPtFuncMC,"MC range MC","l");
    } 
    else h[0]->Draw("e0");
    leg->AddEntry(h[0],"MC Generated tracks ","lep");
    
    if (fPtFunc) 
    {
      fPtFunc->SetLineColor(4); // Blue
      h[0]->Fit(fPtFunc, "IWLR");
      leg->AddEntry(fPtFunc,"MC range Data","l");
    }
  }
  
  if (hRef[0]) // AccEff corrected Data pt histo. 
  {
    hRef[0]->SetLineColor(2); // Red
    if (fPtFuncNew) 
    {
      fPtFuncNew->SetLineColor(2);// Red
      hRef[0]->Fit(fPtFuncNew, "IWLR", "e0sames");
      leg->AddEntry(fPtFuncNew,"Data","l");
    } 
    else hRef[0]->Draw("e0sames");
    leg->AddEntry(hRef[0],"AccEff corrected Data","lep");
  }
  leg->Draw("same");

//------Pad 2
  fcRes->cd(2);
  TLegend*leg2 = new TLegend(0.48,0.7,0.70,0.9);

  if (hAccEff[1]) 
  {
    if (fYFuncMC) 
    {
      fYFuncMC->SetLineColor(3);// Green
      fYFuncMC->SetLineWidth(3);
      h[1]->Fit(fYFuncMC, "IWLR", "e0sames");//  Gen. y histo.
      leg2->AddEntry(fYFuncMC,"MC fit function","l");
    } 
    else h[1]->Draw("");
    leg2->AddEntry(h[1],"MC Generated tracks ","lep");
    
    if (fYFunc) 
    {
      fYFunc->SetLineColor(4);// Blue
      h[1]->Fit(fYFunc, "IWLR");
      leg2->AddEntry(fYFunc,"Old Data fit function","l");
    }
  }
  if (hRef[1])  // AccEff corrected Data y histo. 
  {
    hRef[1]->SetLineColor(2);// Red
    if (fYFuncNew) 
    {
      fYFuncNew->SetLineColor(2); // Red
      hRef[1]->Fit(fYFuncNew, "IWLR", "e0sames");
      leg2->AddEntry(fYFuncNew,"New Data Fit function","l");
    } 
    else hRef[1]->Draw("e0sames");
    leg2->AddEntry(hRef[1],"AccEff corrected Data","lep");
  }
  leg2->Draw("same");

//------Pad 3 and 4
  for (Int_t i = 2; i < 4; i++) 
  {
    fcRes->cd(i+1);
    TLegend*legend = new TLegend(0.48,0.7,0.70,0.9);

    if (i == 2) gPad->SetLogy();
    h[i]->Draw(""); // AccEff corrected Rec.  histo.
    legend->AddEntry(h[i],"Rec. MC Track","lep");
    
    if (hRef[i]) 
    {
      hRef[i]->SetLineColor(2);// Red
      hRef[i]->Draw("sames");// data histo.
      legend->AddEntry(hRef[i],"Data Track","lep");
    }
    legend->Draw("e0sames");
  }
  //__________
  
  printf("Normalizing function....\n");
  
  //__________normalize functions to their integral in the range used in MC 
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
   //__________
  
  // prepare data/MC function ratios
  TF1 *ptRat = (hRat[0] && fPtFunc && fPtFuncNew) ? new TF1("ptRat", this, &AliAnalysisTaskGenTunerJpsi::PtRat, fitRangeMC[0][0], hRat[0]->GetXaxis()->GetXmax(), 0, "AliAnalysisTaskGenTunerJpsi", "PtRat") : 0x0;
  TF1 *yRat = (hRat[1] && fYFunc && fYFuncNew) ? new TF1("yRat", this, &AliAnalysisTaskGenTunerJpsi::YRat, fitRangeMC[1][0], fitRangeMC[1][1], 0, "AliAnalysisTaskGenTunerJpsi", "YRat") : 0x0;
  
  //__________plot ratios
  fcRat = new TCanvas("cRat", "ratios", 900, 600);
  fcRat->Divide(2,2);
  
  for (Int_t i = 0; i < 4 && hRat[i]; i++) 
  {
    fcRat->cd(i+1);
    TLegend*leg= new TLegend(0.48,0.8,0.70,0.9);
    hRat[i]->Draw("e0");
    if (i == 0 && ptRat) 
    {
      printf("ptRat = %p\n", ptRat);
      ptRat->DrawCopy("same");
      leg->AddEntry(ptRat,"MC-Data fit function","l");
      leg->AddEntry(hRat[i],"Data Dist. / MC Gen. ","leg"); 
    }
    else if (i == 1 && yRat)
    {
      yRat->Draw("same");
      leg->AddEntry(ptRat,"MC-Data fit function","l");
      leg->AddEntry(hRat[i],"Data Dist. / MC Gen. ","leg");
    }
    leg->AddEntry(hRat[i],"Data/Rec","leg"); 
    
    leg->Draw("same");
  }
  // __________
  
  //__________print fitting ranges
  if (fPtFuncMC && fYFuncMC) {
    printf("\npT fitting range MC = [%g, %g]\n", fitRangeMC[0][0], fitRangeMC[0][1]);
    printf("y fitting range MC = [%g, %g]\n\n", fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (fPtFunc && fYFunc) {
    printf("pT fitting range = [%g, %g]\n", fitRange[0][0], fitRange[0][1]);
    printf("y fitting range = [%g, %g]\n\n", fitRange[1][0], fitRange[1][1]);
  }
  //__________
  
  //__________print parameters
  if (fPtFuncMC) {
    printf("Double_t oldPtParamMC[%d] = {", fPtFuncMC->GetNpar());
    for (Int_t i = 0; i < fPtFuncMC->GetNpar()-1; i++) printf("%g, ", fPtFuncMC->GetParameter(i));
    printf("%g};\n", fPtFuncMC->GetParameter(fPtFuncMC->GetNpar()-1));
  }
  if (fYFuncMC) {
    printf("Double_t oldYParamMC[%d] = {", fYFuncMC->GetNpar());
    for (Int_t i = 0; i < fYFuncMC->GetNpar()-1; i++) printf("%g, ", fYFuncMC->GetParameter(i));
    printf("%g};\n\n", fYFuncMC->GetParameter(fYFuncMC->GetNpar()-1));
  }
  if (fPtFunc) {
    printf("Double_t oldPtParam[%d] = {", fPtFunc->GetNpar());
    for (Int_t i = 0; i < fPtFunc->GetNpar()-1; i++) printf("%g, ", fPtFunc->GetParameter(i));
    printf("%g};\n", fPtFunc->GetParameter(fPtFunc->GetNpar()-1));
  }
  if (fYFunc) {
    printf("Double_t oldYParam[%d] = {", fYFunc->GetNpar());
    for (Int_t i = 0; i < fYFunc->GetNpar()-1; i++) printf("%g, ", fYFunc->GetParameter(i));
    printf("%g};\n\n", fYFunc->GetParameter(fYFunc->GetNpar()-1));
  }
  if (fPtFuncNew) {
    printf("Double_t newPtParam[%d] = {", fPtFuncNew->GetNpar());
    for (Int_t i = 0; i < fPtFuncNew->GetNpar()-1; i++) printf("%g, ", fPtFuncNew->GetParameter(i));
    printf("%g};\n", fPtFuncNew->GetParameter(fPtFuncNew->GetNpar()-1));
  }
  if (fYFuncNew) {
    printf("Double_t newYParam[%d] = {", fYFuncNew->GetNpar());
    for (Int_t i = 0; i < fYFuncNew->GetNpar()-1; i++) printf("%g, ", fYFuncNew->GetParameter(i));
    printf("%g};\n\n", fYFuncNew->GetParameter(fYFuncNew->GetNpar()-1));
  }
  //__________
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::SetOriginPtFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                              Double_t xMin, Double_t xMax)
{
  /// Create the original function with the parameters used in simulation to generate the pT distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution for cross-check.
  
  assert(param);
  
  delete fPtFuncOld;
  fPtFuncOld = new TF1("fPtFuncOld", formula.Data(), xMin, xMax);
  
  fPtFuncOld->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fPtFuncOld->GetNpar(); ++i)
    if (fixParam[i]) fPtFuncOld->FixParameter(i, fPtFuncOld->GetParameter(i));
  
  NormFunc(fPtFuncOld, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::SetNewPtFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                           Double_t xMin, Double_t xMax)
{
  /// Create the new function with its initial parameters to fit the generated/weighted pT distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution.
  
  assert(param);
  
  delete fPtFuncNew;
  fPtFuncNew = new TF1("fPtFuncNew", formula.Data(), xMin, xMax);
  
  fPtFuncNew->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fPtFuncNew->GetNpar(); ++i)
    if (fixParam[i]) fPtFuncNew->FixParameter(i, fPtFuncNew->GetParameter(i));
  
  NormFunc(fPtFuncNew, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::SetOriginYFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                             Double_t xMin, Double_t xMax)
{
  /// Create the original function with the parameters used in simulation to generate the y distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution for cross-check.
  
  assert(param);
  
  delete fYFuncOld;
  fYFuncOld = new TF1("fYFuncOld", formula.Data(), xMin, xMax);
  
  fYFuncOld->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fYFuncOld->GetNpar(); ++i)
    if (fixParam[i]) fYFuncOld->FixParameter(i, fYFuncOld->GetParameter(i));
  
  NormFunc(fYFuncOld, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::SetNewYFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                          Double_t xMin, Double_t xMax)
{
  /// Create the new function with its initial parameters to fit the generated/weighted y distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution.
  
  assert(param);
  
  delete fYFuncNew;
  fYFuncNew = new TF1("fYFuncNew", formula.Data(), xMin, xMax);
  
  fYFuncNew->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fYFuncNew->GetNpar(); ++i)
    if (fixParam[i]) fYFuncNew->FixParameter(i, fYFuncNew->GetParameter(i));
  
  NormFunc(fYFuncNew, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::SetPtBin(Int_t nofbin,Double_t* bin)
{
  fPtNofBin = nofbin ;

  if(nofbin==0)
  {

    fPtBin = 0x0;
    return;
  } 

  fPtBin = new Double_t[fPtNofBin];

  for (int i = 0; i < fPtNofBin; ++i)
  {
    fPtBin[i]=bin[i];
    // cout << "fPtBin " << i << "=" << fPtBin[i] << endl;
  }
  return;
} 

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::SetYBin(Int_t nofbin,Double_t* bin)
{
  fYNofBin = nofbin;

  if(nofbin==0)
  {
    fYBin = 0x0;
    return;
  } 

  fYBin = new Double_t[fYNofBin];

  for (int i = 0; i < fYNofBin; ++i)
  {
    fYBin[i]=bin[i];
    // cout << "fYBin " << i << "=" << fYBin[i] << endl;
  }
  return;
} 

//________________________________________________________________________
TH1* AliAnalysisTaskGenTunerJpsi::ComputeAccEff(TH1 &hGen, TH1 &hRec, const Char_t *name, const Char_t *title)
{
  /// Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  /// Result is identical to divide histograms with option "B", except here error is forced > 1/gen
  
  Int_t nbins = hGen.GetNbinsX();
  TH1* hAcc = static_cast<TH1*>(hGen.Clone());
  hAcc->SetName(name);
  hAcc->SetTitle(title);
  for (Int_t i = 1; i <nbins+1; i++)
  { 
    Double_t accEff = 0.;
    Double_t accEffErr = 0.;
    Double_t gen = hGen.GetBinContent(i);
    // cout << "gen for bin " << i <<" = " << gen << endl;

    if (gen > 0.) 
    {
      Double_t rec = hRec.GetBinContent(i);
      // cout << "rec for bin " << i <<" = " << rec << endl;

      Double_t genErr = hGen.GetBinError(i);
      // cout << "genErr for bin " << i <<" = " << genErr << endl;

      Double_t recErr = hRec.GetBinError(i);
      // cout << "recErr for bin " << i <<" = " << recErr << endl;


      accEff = rec/gen;
      //Double_t accEffErr2 = ((1.-2.*accEff)*recErr*recErr + accEff*accEff*genErr*genErr)/(gen*gen);
      //accEffErr = TMath::Max(accEff*genErr/gen, TMath::Sqrt(TMath::Abs(accEffErr2)));
      accEffErr = TMath::Max(recErr/gen, accEff*genErr/gen);
      // cout << "AccEff for bin " << i <<" = " << accEff << "+/-" << accEffErr << endl;
    }
    hAcc->SetBinContent(i, accEff);
    hAcc->SetBinError(i, accEffErr);
  }
  
  return hAcc;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTunerJpsi::GetFitLowEdge(TH1 &h)
{
  /// adjust the lower edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindFirstBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0+1)) binAbove0++;
  return h.GetBinLowEdge(binAbove0);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTunerJpsi::GetFitUpEdge(TH1 &h)
{
  /// adjust the upper edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindLastBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0-1)) binAbove0--;
  return h.GetBinLowEdge(binAbove0+1);
}

//________________________________________________________________________
void AliAnalysisTaskGenTunerJpsi::NormFunc(TF1 *f, Double_t min, Double_t max)
{
  /// normalize the function to its integral in the given range
  Double_t integral = f->Integral(min, max);
  if (integral != 0.) f->SetParameter(0, f->GetParameter(0)/integral);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTunerJpsi::PtRat(const Double_t *x, const Double_t */*p*/)
{
  /// generated pT fit function ratio
  // if (fPtFunc ) printf("fPtFunc = %p\n", fPtFunc);
  // if (fPtFuncNew ) printf("fPtFuncNew = %p\n", fPtFuncNew);

  return (fPtFunc && fPtFuncNew) ? fPtFuncNew->Eval(*x) / fPtFunc->Eval(*x) : 0.;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTunerJpsi::YRat(const Double_t *x, const Double_t */*p*/)
{
  /// generated y fit function ratio
  return (fYFunc && fYFuncNew) ? fYFuncNew->Eval(*x) / fYFunc->Eval(*x) : 0.;
}

