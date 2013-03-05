
import FWCore.ParameterSet.Config as cms
from TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff import *


ttbarPhotonMerger = cms.EDFilter("TTGammaMerger",
    ptCut = cms.double(8.),
    drCut = cms.double(.1),
    is2to7 = cms.untracked.bool(True),
    filter = cms.bool(True),
)

ttgammaMerging = cms.Sequence(makeGenEvt * ttbarPhotonMerger)

