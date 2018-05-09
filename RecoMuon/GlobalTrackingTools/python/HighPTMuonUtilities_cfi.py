import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonUpdatorAtVertex_cff import *
from TrackingTools.KalmanUpdators.KFUpdatorESProducer_cfi import *
from TrackingTools.GeomPropagators.SmartPropagator_cff import *
import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi
Chi2EstimatorForMuonTrackLoader = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone()
Chi2EstimatorForMuonTrackLoader.ComponentName = cms.string('Chi2EstimatorForMuonTrackLoader')
Chi2EstimatorForMuonTrackLoader.nSigma = 3.0
Chi2EstimatorForMuonTrackLoader.MaxChi2 = 100000.0

HighPTMuonUtilities = cms.PSet(
		MuonUpdatorAtVertex,


		beamSpot = cms.InputTag('offlineBeamSpot'),

		#Propagator = cms.string('SmartPropagatorRKOpposite'),
		#Propagator = cms.string('SmartPropagatorAnyRK'),
		#Propagator = cms.string('SmartPropagatorAnyOpposite'),
		Propagator = cms.string('SteppingHelixPropagatorOpposite'),

		Selector = cms.string('dxy'),
		#Selector = cms.string('trackRank'),
		trackRankFactor = cms.double(1.0),
		curvPullCut = cms.double(1.)
		#Selector = cms.string('curvPull'),
)

