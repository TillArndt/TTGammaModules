
import FWCore.ParameterSet.Config as cms

ttbarPhotonMerger2to5 = cms.EDFilter("TTGammaMerger",
    ptCut = cms.double(20.),
    drCut = cms.double(.1),
    is2to5 = cms.untracked.bool(True),
    filter = cms.bool(True),
)



