// -*- C++ -*-
//
// Package:    MergingWeightProducer/MergingWeightAnalyzer
// Class:      MergingWeightAnalyzer
// 
/**\class MergingWeightAnalyzer MergingWeightAnalyzer.cc MergingWeightProducer/MergingWeightAnalyzer/plugins/MergingWeightAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrej Saibel
//         Created:  Fri, 24 Feb 2017 14:25:47 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TTree.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MergingWeightAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MergingWeightAnalyzer(const edm::ParameterSet&);
      ~MergingWeightAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
      // Input tags
      const edm::EDGetTokenT<float> TTPlusBMergingWeightToken_;
      
      // Variables in trees
      
      float TTPlusBMergingWeightInTree_;
      
      //Define Tree to fill
      
      TTree*  tree_;
      
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MergingWeightAnalyzer::MergingWeightAnalyzer(const edm::ParameterSet& iConfig): TTPlusBMergingWeightToken_(consumes<float>(iConfig.getParameter<edm::InputTag>("TTPlusBMergingWeight")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


MergingWeightAnalyzer::~MergingWeightAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MergingWeightAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //Get the weight from the miniaod
   edm::Handle<float> TTPlusBMergingWeight;
   iEvent.getByToken(TTPlusBMergingWeightToken_, TTPlusBMergingWeight);
   
   TTPlusBMergingWeightInTree_ = *TTPlusBMergingWeight;
   //fill information into root tree
   
   tree_->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
MergingWeightAnalyzer::beginJob()
{
    edm::Service<TFileService> fileService;
    if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
    
    tree_ = fileService->make<TTree>("tree", "tree");
    tree_->Branch("TTPlusBMergingWeight", &TTPlusBMergingWeightInTree_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MergingWeightAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MergingWeightAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MergingWeightAnalyzer);
