import FWCore.ParameterSet.Config as cms

HighPTMuonRefitter = cms.PSet(


    DTRecSegmentLabel = cms.InputTag("dt1DRecHits"),
    CSCRecSegmentLabel = cms.InputTag("csc2DRecHits"),
    GEMRecHitLabel = cms.InputTag("gemRecHits"),
    ME0RecHitLabel = cms.InputTag("me0Segments"),
    RPCRecSegmentLabel = cms.InputTag("rpcRecHits"),

    MuonHitsOption = cms.int32(1),
    PtCut = cms.double(1.0),
    Chi2ProbabilityCut = cms.double(30.0),
    Chi2CutCSC = cms.double(1.0),
    Chi2CutDT = cms.double(30.0),
    Chi2CutGEM = cms.double(1.0),
    Chi2CutME0 = cms.double(1.0),
    Chi2CutRPC = cms.double(1.0),
    HitThreshold = cms.int32(1),

    Fitter = cms.string('KFFitterForRefitOutsideIn'),
    Smoother = cms.string('KFSmootherForRefitInsideOut'),
    Propagator = cms.string('SmartPropagatorAnyRK'),
    TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
    DoPredictionsOnly = cms.bool(False),
    RefitDirection = cms.string('outsideIn'),
    PropDirForCosmics = cms.bool(False),
    RefitRPCHits = cms.bool(False),
 
    # DYT stuff
    DYTthrs = cms.vint32(10, 10),
    DYTselector = cms.int32(1),
    DYTupdator = cms.bool(True),
    DYTuseAPE = cms.bool(False),


		printStuff = cms.bool( False ),
		minNumHits = cms.int32(1),

)

# This customization will be removed once we get the templates for
# phase2 pixel
from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
phase2_tracker.toModify(HighPTMuonRefitter, TrackerRecHitBuilder = 'WithTrackAngle') # FIXME

