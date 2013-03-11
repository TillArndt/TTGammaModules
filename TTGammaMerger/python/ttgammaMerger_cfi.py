
import FWCore.ParameterSet.Config as cms
from TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff import *


ttbarPhotonMerger2to7 = cms.EDFilter("TTGammaMerger",
    ptCut = cms.double(8.),
    drCut = cms.double(.1),
    is2to7 = cms.untracked.bool(True),
    filter = cms.bool(True),
)
ttgammaMerging2to7 = cms.Sequence(makeGenEvt * ttbarPhotonMerger2to7)

ttbarPhotonMerger2to5 = cms.EDFilter("TTGammaMerger",
    ptCut = cms.double(20.),
    drCut = cms.double(.1),
    is2to5 = cms.untracked.bool(True),
    filter = cms.bool(True),
)
ttgammaMerging2to5 = cms.Sequence(makeGenEvt * ttbarPhotonMerger2to5)



