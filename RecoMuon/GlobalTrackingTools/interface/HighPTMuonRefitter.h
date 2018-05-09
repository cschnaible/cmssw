#ifndef RecoMuon_GlobalTrackingTools_HighPTMuonRefitter_H
#define RecoMuon_GlobalTrackingTools_HighPTMuonRefitter_H

/** \class HighPTMuonRefitter
 *  class to build high pt muon trajectories from 
 *  muon-only information
 *
 *
 *  \author Christian Schnaible (UCLA)
 *
 *
 *
 *
 */

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoMuon/TrackingTools/interface/MuonTrajectoryBuilder.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonReco/interface/DYTInfo.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryParameters.h"
#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonUtilities.h"

namespace edm {class Event;}
namespace reco {class TransientTrack;}

class TrajectoryStateOnSurface;
class TrackerTopology;

class MuonDetLayerMeasurements;
class MuonServiceProxy;

class Trajectory;
class TrajectoryFitter;
class TrajectorySmoother;

class HighPTMuonUtilities;

class HighPTMuonRefitter {

  public:

    typedef TransientTrackingRecHit::RecHitContainer RecHitContainer;
    typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
    typedef TransientTrackingRecHit::RecHitPointer RecHitPointer;
    typedef TransientTrackingRecHit::ConstRecHitPointer ConstRecHitPointer;

    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    typedef std::vector<Trajectory> TC;
    typedef TC::const_iterator TI;

    enum subDetector { PXB = 1, PXF = 2, TIB = 3, TID = 4, TOB = 5, TEC = 6 };

  public:

    /// constructor with Parameter Set and MuonServiceProxy
    HighPTMuonRefitter(const edm::ParameterSet&, const MuonServiceProxy*, const HighPTMuonUtilities*, edm::ConsumesCollector&);
          
    /// destructor
    virtual ~HighPTMuonRefitter();

    /// pass the Event to the algo at each event
    virtual void setEvent(const edm::Event&);

    /// set the services needed by the TrackTransformer
    void setServices(const edm::EventSetup&);

    /// build combined trajectory from subset muon RecHits
    std::pair<Trajectory,Trajectory> refit(
					const reco::Track& globalTrack,
				  const std::string& theMuonHitsOption,
				  const TrackerTopology *tTopo,
					const int &nHits,
				  const int &hitPattern) const;

    /// refit the track with a new set of RecHits
    std::vector<Trajectory> transform(
					const reco::Track& newTrack,
					const reco::TransientTrack track,
					const TransientTrackingRecHit::ConstRecHitContainer& recHitsForReFit) const;
    
    // return DYT-related informations           
    const reco::DYTInfo* getDYTInfo() {return dytInfo;}
    
  protected:

    enum RefitDirection{insideOut,outsideIn,undetermined};
		const MuonServiceProxy* service() const { return theService; }

		// Make a vector<Trajectory> from a free trajectory state and a rechit list
		// Uses propagation only - it assumes the FTS is the correct state
		std::vector<Trajectory> trajFromFTS(
					const Trajectory& traj,
					const FreeTrajectoryState& fts, 
					const ConstRecHitContainer& recHits) const;

		RefitDirection checkRecHitsOrdering(const ConstRecHitContainer&) const;

		// Get the RecHits
		ConstRecHitContainer getRecHits(const reco::TransientTrack& track) const;
		// get rid of Tracker RecHits
		ConstRecHitContainer removeTrackerHits(const ConstRecHitContainer& hits) const;
		// get rid of Muon RecHits
		ConstRecHitContainer removeMuonHits(const ConstRecHitContainer& hits) const;
		// Apply various refit algorithms
		ConstRecHitContainer applyMuonHitsOption(
					const reco::Track& globalTrack, const reco::TransientTrack& track,
					const ConstRecHitContainer& muonRecHits, 
					const std::string& muonHitsOption, const int& nHits, const int& hitpattern,
					std::map<DetId, int> &) const;
			
		// Picky
		/// check muon RecHits, calculate chamber occupancy and select hits to be used in the final fit
		void checkMuonHits(
					const reco::Track&, 
					ConstRecHitContainer&, 
					std::map<DetId, int> &) const;
		/// select muon hits compatible with trajectory; check hits in chambers with showers
		ConstRecHitContainer selectMuonHits(
					const Trajectory&, 
					const std::map<DetId, int> &) const;

		// TPFMS
		/// get the RecHits in the tracker and the first muon chamber with hits 
		ConstRecHitContainer getFirstHits(const ConstRecHitContainer&) const;

		// Combinatoric
		// Select muon hits according to combinatoric pattern
		ConstRecHitContainer applyHitMask(const ConstRecHitContainer&, const int&, const int&) const;
		// Convert hit mask integer to binary
		std::vector<int> convertToBinary(const int&, const int&) const;

		/// Compute chi2 of a Trajectory 
		//double computeChi2(const Trajectory&) const;

		/// print all RecHits of a trajectory
		void printHits(const ConstRecHitContainer&) const;
		void printHitscout(const ConstRecHitContainer&) const;
		void printTrajectories(const Trajectory&, const bool &) const;
		void printTSOS(const TrajectoryStateOnSurface&) const;
		void printRecHit(const TrackingRecHit::ConstRecHitPointer &) const;

  protected:

    std::string theCategory;
    bool theTkTrajsAvailableFlag;
    float thePtCut;

  private:

    const MuonServiceProxy* theService;
		const HighPTMuonUtilities* theHighPTUtilities;
    const edm::Event* theEvent;
  
		std::string   theMuonHitsOption;
    RefitDirection theRefitDirection;

    float theProbCut;
    int   theHitThreshold;
    float theDTChi2Cut;
    float theCSCChi2Cut;
    float theRPCChi2Cut;
    float theGEMChi2Cut;
    float theME0Chi2Cut;
    bool  theCosmicFlag;

    edm::InputTag theME0RecHitLabel;
    edm::Handle<ME0SegmentCollection> theME0RecHits;
    edm::EDGetTokenT<ME0SegmentCollection> theME0RecHitToken;

    edm::InputTag theGEMRecHitLabel;
    edm::Handle<GEMRecHitCollection> theGEMRecHits;
    edm::EDGetTokenT<GEMRecHitCollection> theGEMRecHitToken;

    edm::InputTag theCSCRecHitLabel;
    edm::Handle<CSCRecHit2DCollection> theCSCRecHits;
    edm::EDGetTokenT<CSCRecHit2DCollection> theCSCRecHitToken;

    edm::InputTag theDTRecHitLabel;
    edm::EDGetTokenT<DTRecHitCollection> theDTRecHitToken;
    edm::Handle<DTRecHitCollection>    theDTRecHits;

    edm::Handle<CSCSegmentCollection> CSCSegments;
    edm::EDGetTokenT<CSCSegmentCollection> CSCSegmentsToken;

    edm::Handle<DTRecSegment4DCollection> all4DSegments;
    edm::EDGetTokenT<DTRecSegment4DCollection> all4DSegmentsToken;

    unsigned long long theCacheId_TRH;        

    bool theRPCInTheFit;
		bool printStuff;
		unsigned int minNumHits;
    double theRescaleErrorFactor;

    std::vector<int> theDYTthrs;
    int theDYTselector;
    bool theDYTupdator;
    bool theDYTuseAPE;
    reco::DYTInfo *dytInfo;

    std::string thePropagatorName;
    std::string theFitterName;
    std::unique_ptr<TrajectoryFitter> theFitter;
    std::string theSmootherName;
    std::unique_ptr<TrajectorySmoother> theSmoother;
  
    std::string theTrackerRecHitBuilderName;
    edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
    TkClonerImpl hitCloner;
  
    std::string theMuonRecHitBuilderName;
    edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;


};
#endif
