#ifndef RecoMuon_GlobalTrackingTools_HighPTMuonUtilities_H
#define RecoMuon_GlobalTrackingTools_HighPTMuonUtilities_H

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToValue.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

namespace edm {class Event; class EventSetup; class ParameterSet;}

class Trajectory;
class Propagator;
class TrajectorySmoother;
class MuonServiceProxy;
class MuonUpdatorAtVertex;

class HighPTMuonUtilities { 
	public:
	typedef ROOT::Math::SMatrix< double,10,10,ROOT::Math::MatRepSym<double,  10> > Matrix1010;
	typedef ROOT::Math::SMatrix< double,10, 5,ROOT::Math::MatRepStd<double,10,5> > Matrix105;
	typedef ROOT::Math::SMatrix< double, 5,10,ROOT::Math::MatRepStd<double,5,10> > Matrix510;
	typedef ROOT::Math::SMatrix< double, 5, 5,ROOT::Math::MatRepStd<double,5, 5> > Matrix55;
	typedef ROOT::Math::SVector< double,10> Vector10;

    /// constructor with Parameter Set and MuonServiceProxy
    HighPTMuonUtilities(const edm::ParameterSet&, const MuonServiceProxy*, edm::ConsumesCollector&);
          
    /// destructor
    virtual ~HighPTMuonUtilities();
		//
    /// pass the Event to the algo at each event
    virtual void setEvent(const edm::Event&);

		FreeTrajectoryState updateMuonOnlyTraj(
				const Trajectory& refitTraj,
				const reco::Track& globalTrack) const;

		std::pair<Trajectory,Trajectory> select(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack,
			const reco::GenParticle& genMuon) const;
			//const reco::Track& glbTrack) const;

		reco::TrackExtra buildTrackExtra(const Trajectory& trajectory) const;
		std::pair<bool,reco::Track> buildTrackFromTrajAtPCA(
				const Trajectory&) const;
		std::pair<bool,reco::Track> buildTrackFromTrajAtPCA(
				const Trajectory&, const reco::BeamSpot&, const bool& vtxUpdate) const;
		std::pair<bool,reco::Track> buildTrackFromTrackAtPCA(
				const reco::Track&, const reco::BeamSpot&) const;

	private:

    const edm::Event* theEvent;
    const MuonServiceProxy *theService;
		MuonUpdatorAtVertex *theUpdatorAtVtx;

		std::string thePropagatorName;
		std::string theSelectorName;
		double factor;
		double curvPullCut;

		edm::EDGetTokenT<reco::BeamSpot> theBeamSpotToken;
		edm::InputTag theBeamSpotInputTag;
		edm::Handle<reco::BeamSpot> beamSpot;

		int calculateNDoF(const Trajectory& traj, const bool& bon, const bool& vtxUpdate) const;

		FreeTrajectoryState buildFTSfromParCov(
				const std::pair<CurvilinearTrajectoryParameters,CurvilinearTrajectoryError>& refitUpdateParCov,
				const reco::Track& globalTrackAtPCA) const;

    std::pair<CurvilinearTrajectoryParameters,CurvilinearTrajectoryError> KFupdateTrackWithVtx(
				const reco::Track &refit, const reco::Track &glbTrack) const;

		std::pair<Trajectory,Trajectory> selectBasedOnDxyPull(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack) const;
			//const reco::Track& glbTrack) const;

		std::pair<Trajectory,Trajectory> selectBasedOnTrackRank(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack) const;
			//const reco::Track& glbTrack) const;

		std::pair<Trajectory,Trajectory> selectBasedOnCurvPull(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack,
			const std::string& refTrack) const;
			//const reco::Track& glbTrack) const;

		std::pair<Trajectory,Trajectory> selectBasedOnRelCurvErr(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack) const;

		std::pair<Trajectory,Trajectory> selectBasedOnGEN(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack,
			const reco::GenParticle& genMuon) const;
		std::pair<Trajectory,Trajectory> selectBasedOnTEST(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack,
			const reco::GenParticle& genMuon) const;
			//const reco::Track& glbTrack) const;
		std::pair<Trajectory,Trajectory> selectBasedOnTEST_2(
			const std::vector< std::pair<Trajectory,Trajectory> >& refits,
			const reco::Muon& glbTrack,
			const reco::GenParticle& genMuon) const;
			//const reco::Track& glbTrack) const;

		void combine(const reco::Track &tracker, const reco::Track &refit) const;
		std::pair<double,double> combineRet(const reco::Track &tracker, const reco::Track &refit) const;


		void printTSOS(const TrajectoryStateOnSurface& thisTSOS) const;
		template<class T> void printTrack(const T& track, const std::string& name) const;
		template<class T> void printComb(const T& track, const Trajectory& traj, const std::string& name, const int& nHits, const int& hitOption) const;
		std::string hitMaskPrint(const int& nHits, const int& hitOption) const;

};
#endif
