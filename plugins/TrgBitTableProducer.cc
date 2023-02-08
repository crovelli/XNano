//// table to produce hlt bits that we are going to use in the analysis, to avoid saving every bit in the menu

// system include files
#include <memory>

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TString.h"
#include <string>

class TrgBitTableProducer : public edm::stream::EDProducer<> {

public:

  explicit TrgBitTableProducer(const edm::ParameterSet &cfg):
    hltresultsToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag> ("hltresults"))),
    hltpaths_( cfg.getParameter< std::vector<std::string> >( "paths" ) )
  {
    produces<nanoaod::FlatTable>();
  }
  
  ~TrgBitTableProducer() override {}
  
  void produce(edm::Event&, edm::EventSetup const&) override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  
  const edm::EDGetTokenT< edm::TriggerResults >    hltresultsToken_;
  const std::vector< std::string >                 hltpaths_;
  TString * algoBitToName  =  new TString[512]; 
  bool loaded = false;
};

void 
TrgBitTableProducer::produce( edm::Event &evt, edm::EventSetup const &stp) 
{
  edm::Handle< edm::TriggerResults > hltResults;
  evt.getByToken( hltresultsToken_, hltResults);
  
  std::vector<int> hltbits;
  unsigned int Npaths = hltpaths_.size();
  hltbits.reserve( Npaths );
  edm::TriggerNames trigName = evt.triggerNames( *hltResults );   
  
  // get HLT triggers
  if ( hltResults.failedToGet() ){
    for ( unsigned int ibit = 0; ibit < Npaths; ++ibit)
      hltbits.push_back( 0 );
    
  } else {
    int Ntrg = hltResults->size();
    for ( auto& hltpath: hltpaths_ ){
      bool fire = false; 
      for( int itrg = 0; itrg < Ntrg; ++itrg ){
        if ( !hltResults->accept( itrg ) ) continue;
        TString TrigPath = trigName.triggerName( itrg );
        if ( TrigPath.Contains( hltpath ) ){
	  fire=true;
	  break; 
        }
      }
    
      if( fire ) hltbits.push_back( 1 );
      else hltbits.push_back( 0 );
    }
  }
  
  auto tab  = std::make_unique<nanoaod::FlatTable>(1,"", true);
  for (unsigned int ipath = 0; ipath <Npaths; ++ipath ){
    //std::cout << "ipath = " << ipath << ", hltpaths_[ipath] = " << hltpaths_[ipath] << ", hltbits[ipath] = " << hltbits[ipath] << std::endl; 
    tab->addColumnValue<int> (hltpaths_[ipath], hltbits[ipath], "hlt path");
  }

  evt.put(std::move(tab));
}

// define this as a plug-in
DEFINE_FWK_MODULE(TrgBitTableProducer);
