#ifndef ALIANALYSISMUMUMIXEVT_H
#define ALIANALYSISMUMUMIXEVT_H

/**
 *
 * \class AliAnalysisMuMuNch
 * \brief Invariant mass dimuon analysis
 * \author L. Aphecetche and J. Martin Blanco (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "TString.h"
#include "TH2.h"

class TH2F;
class AliVParticle;

class AliAnalysisMuMuMixEvt : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuMixEvt();
  virtual ~AliAnalysisMuMuMixEvt();
  
  Bool_t IsPtInRange(const AliVParticle& t1, const AliVParticle& t2,
                           Double_t& ptmin, Double_t& ptmax) const;
  
  void NameOfIsPtInRange(TString& name, Double_t& ymin, Double_t& ymax) const;
  
  Bool_t IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2) const;
  void NameOfIsRapidityInRange(TString& name) const { name = "PAIRY"; }
    
protected:
  
  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality);

  virtual void FillHistosMixEventForPair(const char* eventSelection,const char* triggerClassName,
                                 const char* centrality,
                                 const char* pairCutName,
                                 const AliVParticle& part,
                                 const AliVParticle& part2);
    
private:
  
  void CreateMinvHistograms(const char* eventSelection, const char* triggerClassName, const char* centrality);
  
  TString GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const;

private:
  
  Int_t fsystLevel;
  
  ClassDef(AliAnalysisMuMuMixEvt,1) // implementation of AliAnalysisMuMuBase for muon pairs
};

#endif

