import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *
from RecoMuon.GlobalTrackingTools.GlobalTrajectoryBuilderCommon_cff import *
from RecoMuon.GlobalTrackingTools.HighPTMuonRefitter_cff import *
from RecoMuon.GlobalTrackingTools.HighPTMuonUtilities_cfi import *
highPTMuons = cms.EDProducer("HighPTMuonProducer",
    MuonTrackLoaderForGLB,
    #    InputTag MuonCollectionLabel = standAloneMuons:UpdatedAtVtx
    MuonServiceProxy,
    MuonCollectionLabel = cms.InputTag("globalMuons"),
    Refits = cms.vstring(
		'default', 
		'firstHit', 
		'picky', 
		'dyt',
		'combinatoric',
		'tracker'
	),
    RefitterParameters = cms.PSet(
       HighPTMuonRefitter 
    ),
	UtilitiesParameters = cms.PSet(
		HighPTMuonUtilities
	),
)

