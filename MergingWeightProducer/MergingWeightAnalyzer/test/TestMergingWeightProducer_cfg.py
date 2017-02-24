import FWCore.ParameterSet.Config as cms

process = cms.Process("Analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)

process.source = cms.Source("PoolSource",
    # insert any 4 or 5 flavor scheme miniaod file that you want to test
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring15MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/06D46D97-C66D-E511-9ABF-00266CFAE20C.root'
    )
)
    
    
#Definitions that are required for the genHFHadronMatcher
#as copied from genTTbarCategorizer
genParticleCollection = 'prunedGenParticles'
genJetInputParticleCollection = 'slimmedGenJets'
genJetCollection = 'slimmedGenJets'
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")



from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu.clone(
	src = genJetInputParticleCollection
)

## Produce own jets (re-clustering in miniAOD needed at present to avoid crash)
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsCustom = ak4GenJets.clone(
    src = 'genParticlesForJetsNoNu',
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt")
)

## Ghost particle collection used for Hadron-Jet association
# MUST use proper input particle collection
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection,
)
#Input particle collection for matching to gen jets (partons + leptons)
# MUST use use proper input jet collection: the jets to which hadrons should be associated
# rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
# More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition
from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
)

# Plugin for analysing B hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
    genParticles = genParticleCollection,
    jetFlavourInfos = "genJetFlavourInfos",
    onlyJetClusteredHadrons = cms.bool(False)
)

from MergingWeightProducer.MergingWeightProducer.MergingWeightProducer_cfi import MergingWeightProducer
process.MergingWeightProducer = MergingWeightProducer.clone(
        genJets = cms.InputTag(genJetCollection)
    
)

process.MergingWeightAnalyzer = cms.EDAnalyzer('MergingWeightAnalyzer',
                TTPlusBMergingWeight = cms.InputTag("MergingWeightProducer","TTPlusBMergingWeight")
)

#register output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("TTPlusBMergingWeight.root")
)

process.p = cms.Path(process.matchGenBHadron*process.MergingWeightProducer*process.MergingWeightAnalyzer)
