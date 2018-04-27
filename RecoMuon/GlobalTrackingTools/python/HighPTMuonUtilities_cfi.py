import FWCore.ParameterSet.Config as cms

from TrackingTools.KalmanUpdators.KFUpdatorESProducer_cfi import *
from TrackingTools.GeomPropagators.SmartPropagator_cff import *
import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi
Chi2EstimatorForMuonTrackLoader = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone()
Chi2EstimatorForMuonTrackLoader.ComponentName = cms.string('Chi2EstimatorForMuonTrackLoader')
Chi2EstimatorForMuonTrackLoader.nSigma = 3.0
Chi2EstimatorForMuonTrackLoader.MaxChi2 = 100000.0

HighPTMuonUtilities = cms.PSet(

		beamSpot = cms.InputTag('offlineBeamSpot'),

		#Propagator = cms.string('SmartPropagatorAnyRK'),
		Propagator = cms.string('SteppingHelixPropagatorOpposite'),

		Selector = cms.string('dxy'),
		#Selector = cms.string('trackRank'),
		trackRankFactor = cms.double(1.0),
)

