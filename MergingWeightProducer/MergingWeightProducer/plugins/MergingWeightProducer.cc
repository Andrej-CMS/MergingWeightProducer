// -*- C++ -*-
//
// Package:    MergingWeightProducer/MergingWeightProducer
// Class:      MergingWeightProducer
// 
/**\class MergingWeightProducer MergingWeightProducer.cc MergingWeightProducer/MergingWeightProducer/plugins/MergingWeightProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrej Saibel
//         Created:  Thu, 23 Feb 2017 10:19:36 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"


//
// class declaration
//

class MergingWeightProducer : public edm::stream::EDProducer<> {
   public:
      explicit MergingWeightProducer(const edm::ParameterSet&);
      ~MergingWeightProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      //virtual double CalculateMergingWeight(const double, const double, const double);
      //virtual double GetHardestBGenJetPt(const std::vector<reco::GenJetCollection>& , const std::vector<int>& , const std::vector<int>& , const std::vector<int>&);
      
      double PtMin_;
      double PtMax_;
      
      
      //InputTags
      const edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
      
      const edm::EDGetTokenT<std::vector<int> > genBHadIndexToken_;
      const edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken_;
      const edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken_;

      // ----------member data ---------------------------
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
MergingWeightProducer::MergingWeightProducer(const edm::ParameterSet& iConfig): genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
genBHadIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadIndex"))),
genBHadFlavourToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadFlavour"))),
genBHadJetIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadJetIndex")))
{
   //register your products
    produces<float>("TTPlusBMergingWeight");
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
   PtMin_ = iConfig.getParameter<double>("PtMin");
   PtMax_ = iConfig.getParameter<double>("PtMax");
  
}


MergingWeightProducer::~MergingWeightProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MergingWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    //using namespace edm;
    edm::Handle<reco::GenJetCollection> genJets;    
    iEvent.getByToken(genJetsToken_, genJets);
    
    edm::Handle<std::vector<int> > genBHadIndex;
    iEvent.getByToken(genBHadIndexToken_, genBHadIndex);
    
    edm::Handle<std::vector<int> > genBHadFlavour;
    iEvent.getByToken(genBHadFlavourToken_, genBHadFlavour);
    
    edm::Handle<std::vector<int> > genBHadJetIndex;
    iEvent.getByToken(genBHadJetIndexToken_, genBHadJetIndex);
    
    double HardestBGenJetPt = -1.0;
    std::cout << "hello does this compile " << genBHadIndex->size() << std::endl;
    //find additional b hadrons and jets
    for(size_t BHadronID = 0; BHadronID < genBHadIndex->size(); ++BHadronID){
        const int BJetIndex = genBHadJetIndex->at(BHadronID);
        std::cout << "----------------Event---------------" << std::endl;
        //loop only over jets with associated b-hadrons
        if(BJetIndex < 0 ) continue;
        
        const int BHadronFlavour = genBHadFlavour->at(BHadronID);
        std::cout << "HadronID : " << BHadronID << " , BHadronFlavour : " << BHadronFlavour << " , BJetIndex : " << BJetIndex << std::endl;
        //a b hadron/jet is concidered to be additional, if it doesnt come from the decay of the Top, W, Z, or Higgs
        if(std::abs(BHadronFlavour)!=6 && std::abs(BHadronFlavour)!=24 && std::abs(BHadronFlavour)!=23 && std::abs(BHadronFlavour)!=25){
            if(std::abs(genJets->at(BJetIndex).pt()) > HardestBGenJetPt)
                HardestBGenJetPt = std::abs(genJets->at(BJetIndex).pt());
            
        }
        
    }
    
    //const float HardestBGenJetPt = GetHardestBGenJetPt(*genJets, *genBHadIndex, *genBHadFlavour, *genBHadJetIndex);
    
    //const float TTPlusBMergingWeight = CalculateMergingWeight( HardestBGenJetPt, PtMin_, PtMax_) ;
    std::cout << "HardestBGenJetPt : " << HardestBGenJetPt << std::endl;
    float Weight = -1.0;
    
    if(HardestBGenJetPt < PtMin_){
       Weight = 0.0;
    }
    else if((HardestBGenJetPt >= PtMin_) && (HardestBGenJetPt <= PtMax_)) {
        Weight = 0.5*(1.0-std::cos(float(ROOT::Math::Pi())*(HardestBGenJetPt-PtMin_)/(PtMax_+PtMin_)));
    }
    else  Weight = 1.0;
    
    std::cout << "Weight : " << Weight << std::endl;
    
    std::auto_ptr<float> TTPlusBMergingWeight(new float);
    
    *TTPlusBMergingWeight = Weight;
    
    std::cout << "TTPlusBMergingWeight : "  << *TTPlusBMergingWeight << std::endl;
    iEvent.put(TTPlusBMergingWeight,"TTPlusBMergingWeight");
   
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MergingWeightProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MergingWeightProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MergingWeightProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MergingWeightProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MergingWeightProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MergingWeightProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
/*float GetHardestBGenJetPt(const std::vector<reco::GenJet>& genJets, const std::vector<int>& genBHadIndex, const std::vector<int>& genBHadFlavour, const std::vector<int>& genBHadJetIndex){
    float HardestBGenJetPt = -1.0;
    //find additional b hadrons and jets
    for(uint BHadronID = 0; BHadronID < genBHadIndex.size(); BHadronID++){
        const int BJetIndex = genBHadJetIndex[BHadronID];
        //loop only over jets with associated b-hadrons
        if(BJetIndex < 0 ) continue;
        
        const int BHadronFlavour = genBHadFlavour[BHadronID];
        //a b hadron/jet is concidered to be additional, if it doesnt come from the decay of the Top, W, Z, or Higgs
        if(std::abs(BHadronFlavour)!=6 && std::abs(BHadronFlavour)!=24 && std::abs(BHadronFlavour)!=23 && std::abs(BHadronFlavour)!=25){
            if(std::abs(genJets[BJetIndex].pt()) > HardestBGenJetPt)
                HardestBGenJetPt = std::abs(genJets[BJetIndex].pt());
            
        }
    }
    
    return HardestBGenJetPt;
}

float CalculateMergingWeight(const double HardestBGenJetPt, const double PtMin, const double PtMax) {
    if(HardestBGenJetPt < PtMin){
        return 0.0;
    }
    else if((HardestBGenJetPt >= PtMin) && (HardestBGenJetPt <= PtMax)) {
        return 0.5*(1.0-std::cos(float(ROOT::Math::Pi())*(HardestBGenJetPt-PtMin)/(PtMax+PtMin)));
    }
    else 
        return 1.0;
}*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MergingWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MergingWeightProducer);
