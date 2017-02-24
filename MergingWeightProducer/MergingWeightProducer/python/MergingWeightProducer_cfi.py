import FWCore.ParameterSet.Config as cms

MergingWeightProducer = cms.EDProducer('MergingWeightProducer',
            #iputtags from root file (Mini)AOD                           
                                       
            #iConfig parameters
            PtMin = cms.double(15.0),
            PtMax = cms.double(25.0),
            #InputTags
            genJets = cms.InputTag("ak4GenJets"),
            genBHadIndex = cms.InputTag("matchGenBHadron", "genBHadIndex"),
            genBHadFlavour = cms.InputTag("matchGenBHadron", "genBHadFlavour"),
            genBHadJetIndex = cms.InputTag("matchGenBHadron", "genBHadJetIndex")
)
