#include "AliAnalysisMuMuMixEvt.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuMixEvt
 *
 * Analysis which fills a bunch of mixing event histograms for invariant mass analysis of J/psi
 *
 * Can be used on real data and/or simulated (for instance to get Acc x Eff)
 *
 * Can optionally use as input an already computed Acc x Eff matrix that will be applied
 * when filling the invariant mass histograms.
 *
 */

#include "TH2F.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "AliAnalysisMuMuBinning.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "AliMCEvent.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"
#include "TParameter.h"
#include <cassert>

ClassImp(AliAnalysisMuMuMixEvt)

//_____________________________________________________________________________
AliAnalysisMuMuMixEvt::AliAnalysisMuMuMixEvt()
: AliAnalysisMuMuBase(),
{

}

//_____________________________________________________________________________
AliAnalysisMuMuMixEvt::~AliAnalysisMuMuMixEvt()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliAnalysisMuMuMixEvt::DefineHistogramCollection(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality)
{
  /// Define the histograms this analysis will use
  
  // Check if histo is not already here
  if ( Histo(eventSelection,triggerClassName,centrality,"AliAnalysisMuMuMixEvt") )
  {
    return;
  }
  
  // dummy histogram to signal that we already defined all our histograms (see above)
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"AliAnalysisMuMuMixEvt","Dummy semaphore",1,0,1);

  /// Create invariant mass histograms
  
  //__________mass range
  Double_t minvMin = 2;
  Double_t minvMax = 8;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,0.025);
  
  Int_t nMCMinvBins = GetNbins(minvMin,minvMax,0.1);
  //__________
  
  //__________Rapidity range
  Double_t rapidityMin = -5;
  Double_t rapidityMax = -2;
  Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);
  //__________
  
  //__________eta range
  Double_t etaMin = -5;
  Double_t etaMax = -2;
  Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);
  

  //__________Set histo
  
  //___Histos for pair
  TObjArray* bins = Binning()->CreateBinObjArray("psi","integrated,ptvsy,yvspt,pt,y,phi","");//We may include ,v0a,v0acent

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","#mu+#mu- Pt mixed distribution",
                  nMinvBins,minvMax,minvMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","#mu+#mu- Y mixed distribution",
                   nbinsRapidity,rapidityMin,rapidityMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","#mu+#mu- Eta mixed distribution",
                   nbinsEta,etaMin,etaMax); 

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","#mu+#mu+ Pt mixed distribution",
                  nMinvBins,minvMax,minvMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","#mu+#mu+ Y mixed distribution",
                   nbinsRapidity,rapidityMin,rapidityMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","#mu+#mu+ Eta mixed distribution",
                   nbinsEta,etaMin,etaMax);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","#mu-#mu- Pt mixed distribution",
                  nMinvBins,minvMax,minvMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","#mu-#mu- Y mixed distribution",
                   nbinsRapidity,rapidityMin,rapidityMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","#mu-#mu- Eta mixed distribution",
                   nbinsEta,etaMin,etaMax);
  //____
  
  //___Histos for pure MC
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","MCINPUT #mu+#mu- Pt mixed distribution",
                   nMinvBins,minvMax,minvMax,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","MCINPUT #mu+#mu- Y mixed distribution",
                   nbinsRapidity,rapidityMin,rapidityMax,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","MCINPUT #mu+#mu- Eta mixed distribution",
                   nbinsEta,etaMin,etaMax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Pt","MCINPUT #mu+#mu- Pt mixed distribution",
                    nMinvBins,minvMax,minvMax,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Y","MCINPUT #mu+#mu- Y mixed distribution",
                    nbinsRapidity,rapidityMin,rapidityMax,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Eta","MCINPUT #mu+#mu- Eta mixed distribution",
                    nbinsEta,etaMin,etaMax);
  //____

  TIter next(bins);
  AliAnalysisMuMuBinning::Range* r;
  Int_t nb(0);
  //__________Create Minv Histos for each bin 
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    TString minvName(GetMinvHistoName(*r,kFALSE));// Histos name
    
    ++nb;
    
    if ( !IsHistogramDisabled(minvName.Data()) )
    {
      
      AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));
      
      CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                       Form("#mu+#mu- mixed event inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                            r->AsString().Data()),
                       nMinvBins,minvMin,minvMax,-2);
      
      CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                       Form("MCINPUT #mu+#mu- mixed event inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                            r->AsString().Data()),
                       nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo
      
      CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),minvName.Data(),
                        Form("MCINPUT #mu+#mu- mixed event inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                             r->AsString().Data()),
                        nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo 
    }
  }

  delete bins;
}

//_____________________________________________________________________________
void AliAnalysisMuMuMixEvt::FillHistosMixEventForPair(const char* eventSelection, const char* triggerClassName,
                                            const char* centrality, const char* pairCutName,
                                            const AliVParticle& tracki,
                                            const AliVParticle& trackj)
{
  /// Fill pair Histos for given event / trigger / centrality / pair cut.
  /// Works eather on MC or Data.
  
  // 4-Vector
  TLorentzVector pi(tracki.Px(),tracki.Py(),tracki.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+tracki.P()*tracki.P()));
  // Select only muons
  if (!AliAnalysisMuonUtility::IsMuonTrack(&tracki) ) return;
  if (!AliAnalysisMuonUtility::IsMuonTrack(&trackj) ) return;
  
  TLorentzVector pair4Momentum(trackj.Px(),trackj.Py(),trackj.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+trackj.P()*trackj.P()));
  
  pair4Momentum += pi;

  // Select pair charge
  Float_t PairCharge = tracki.Charge() + trackj.Charge();
  char* sPaireCharge = Form("PaireCharge=%f",PairCharge);

  Double_t inputWeight = 1.;
  
  // ________Fill histo
  if ( !IsHistogramDisabled("Pt") )
  {
    Histo(eventSelection,triggerClassName,centrality,pairCutName,sPaireCharge,"Pt")->Fill(pair4Momentum.Pt(),inputWeight);
  }
  if ( !IsHistogramDisabled("Y") )
  {
    Histo(eventSelection,triggerClassName,centrality,pairCutName,sPaireCharge,"Y")->Fill(pair4Momentum.Rapidity(),inputWeight);
  }
  if ( !IsHistogramDisabled("Eta") )
  {
    Histo(eventSelection,triggerClassName,centrality,pairCutName,sPaireCharge,"Eta")->Fill(pair4Momentum.Eta());
  }
  // ________
  

  
  TObjArray* bins = Binning()->CreateBinObjArray("psi","integrated,ptvsy,yvspt,pt,y,phi",""); // We may include: ,v0a,v0acent
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* r;

  // ----- Loop over each bin type. Check if pair pass all tests and fill histo (except pt an y alone) -----
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    Bool_t ok(kFALSE);

   if ( r->Is2D() ) // Bin x vs Bin y histos tests
    {
      if ( r->AsString().BeginsWith("PTVSY") )
      {
        ok = r->IsInRange(pair4Momentum.Rapidity(),pair4Momentum.Pt());
      }
      else if ( r->AsString().BeginsWith("YVSPT") )
      {
        ok = r->IsInRange(pair4Momentum.Pt(),pair4Momentum.Rapidity());
      }
      else
      {
        AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
      }
    }
    else // the rest
    {
      if ( r->Quantity() == "PT" )
      {
        ok = r->IsInRange(pair4Momentum.Pt());
      }
      else if ( r->Quantity() == "Y" )
      {
        ok = r->IsInRange(pair4Momentum.Rapidity());
      }
      else if ( r->Quantity() == "PHI" )
      {
        ok = r->IsInRange(pair4Momentum.Phi());
      }
    }
    
    if ( ok ) // Pair pass all conditions, either MC or not
    {
      TString minvName = GetMinvHistoName(*r,sPaireCharge);
      
      if (!IsHistogramDisabled(minvName.Data()))
      {
        TH1* h(0x0);
        if ( ok )
        {
          h = Histo(eventSelection,triggerClassName,centrality,pairCutName,minvName.Data());
          
          if (!h)
          {
            AliError(Form("Could not get %s",minvName.Data()));
            //continue;
          }
          else h->Fill(pair4Momentum.M(),inputWeight);// Fill Minv. Histo
        }
      }
    }
  }
  //_________________________
  
  delete bins;
  
  
}


//_____________________________________________________________________________
TString AliAnalysisMuMuMixEvt::GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, char* charge) const
{
  return TString::Format("MinvMixedUS%s%s",r.AsString().Data(),
                         charge);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMixEvt::IsPtInRange(const AliVParticle& t1, const AliVParticle& t2, Double_t& ptmin, Double_t& ptmax) const
{
  /// Whether the pair passes the rapidity cut
  
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t pt = total.Pt();
  
  return  ( pt < ptmax && pt > ptmin );
}

//_____________________________________________________________________________
void AliAnalysisMuMuMixEvt::NameOfIsPtInRange(TString& name, Double_t& ptmin, Double_t& ptmax) const
{
  name.Form("PAIRPTIN%2.1f-%2.1f",ptmin,ptmax);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMixEvt::IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2) const
{
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t y = total.Rapidity();

  return  ( y < -2.5 && y > -4.0 );
}
