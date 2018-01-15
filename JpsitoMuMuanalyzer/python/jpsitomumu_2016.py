print "\n=> running on 2016 data \n"

#####################
#  cmssw configs    #
#####################

import FWCore.ParameterSet.Config as cms
from jpsitomumu_2016_cfi import process 

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1002) )

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(#'/store/data/Run2016F/Charmonium/AOD/23Sep2016-v1/50000/02ABBB71-D495-E611-A08A-0CC47A4D76BE.root')
#'/store/data/Run2016H/Charmonium/AOD/PromptReco-v2/000/281/207/00000/0042DCF9-6382-E611-ABFD-02163E0133B7.root')
#'/store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/273/158/00000/14E579B3-271A-E611-911E-02163E013584.root')
#'/store/data/Run2016B/DoubleMuonLowMass/AOD/23Sep2016-v3/00000/001AD4E7-4498-E611-B144-002590E7DF2A.root')
'/store/data/Run2016G/DoubleMuonLowMass/AOD/23Sep2016-v1/90000/F2D4D1D6-1597-E611-8963-0CC47A0AD63E.root')
)

#process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data') 

# do trigger matching for muons
triggerProcessName = 'HLT'
print "\n line 27"
process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
    # match by DeltaR only (best match by DeltaR)
    'PATTriggerMatcherDRLessByR',                         
    src                   = cms.InputTag('cleanPatMuons'),
# default producer label as defined in
# PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4_3_Jpsi_Displaced*",0,0)'),
    maxDeltaR             = cms.double(0.1),
    # only one match per trigger object
    resolveAmbiguities    = cms.bool(True),
    # take best match found per reco object (by DeltaR here, see above)       
    resolveByMatchQuality = cms.bool(False))

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'],
                              hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
    ('HLT_DoubleMu4_3_Jpsi_Displaced',  'hltDisplacedmumuFilterDoubleMu43Jpsi'),
    ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
