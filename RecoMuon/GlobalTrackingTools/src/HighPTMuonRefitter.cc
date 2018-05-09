/**
 *  Class: HighPTMuonRefitter
 *
 *  Description:
 *  Refit global muon with muon information only
 *  Apply global muon vertex position and direction constraint
 *
 *
 *  \author	Christian Schnaible (UCLA)
 *
 **/

#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonRefitter.h"

//---------------
// C++ Headers --
//---------------

#include <iostream>
#include <iomanip>
#include <algorithm>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include "TrackingTools/TrackFitters/interface/RecHitLessByDet.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/ME0Segment.h>
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "RecoMuon/GlobalTrackingTools/interface/DynamicTruncation.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHitBuilder.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

using namespace std;
using namespace edm;

//----------------
// Constructors --
//----------------

HighPTMuonRefitter::HighPTMuonRefitter(const edm::ParameterSet& par,
				       const MuonServiceProxy* service,
							 const HighPTMuonUtilities* highPTUtilities,
				       edm::ConsumesCollector& iC) : 
  theService(service),
	theHighPTUtilities(highPTUtilities),
  theCosmicFlag(par.getParameter<bool>("PropDirForCosmics")),
  theME0RecHitLabel(par.getParameter<InputTag>("ME0RecHitLabel")),
  theGEMRecHitLabel(par.getParameter<InputTag>("GEMRecHitLabel")),
  theCSCRecHitLabel(par.getParameter<InputTag>("CSCRecSegmentLabel")),
  theDTRecHitLabel(par.getParameter<InputTag>("DTRecSegmentLabel")) {

  theCategory = par.getUntrackedParameter<string>("Category", "Muon|RecoMuon|GlobalMuon|HighPTMuonRefitter");

  theHitThreshold = par.getParameter<int>("HitThreshold");
  theDTChi2Cut  = par.getParameter<double>("Chi2CutDT");
  theCSCChi2Cut = par.getParameter<double>("Chi2CutCSC");
  theRPCChi2Cut = par.getParameter<double>("Chi2CutRPC");
  theGEMChi2Cut = par.getParameter<double>("Chi2CutGEM");
  theME0Chi2Cut = par.getParameter<double>("Chi2CutME0");

  // Refit direction
  string refitDirectionName = par.getParameter<string>("RefitDirection");

  if (refitDirectionName == "insideOut" ) theRefitDirection = insideOut;
  else if (refitDirectionName == "outsideIn" ) theRefitDirection = outsideIn;
  else 
    throw cms::Exception("TrackTransformer constructor") 
      <<"Wrong refit direction chosen in TrackTransformer ParameterSet"
      << "\n"
      << "Possible choices are:"
      << "\n"
      << "RefitDirection = insideOut or RefitDirection = outsideIn";
  
  theFitterName = par.getParameter<string>("Fitter");  
  theSmootherName = par.getParameter<string>("Smoother"),  
  thePropagatorName = par.getParameter<string>("Propagator");


  theTrackerRecHitBuilderName = par.getParameter<string>("TrackerRecHitBuilder");
  theMuonRecHitBuilderName = par.getParameter<string>("MuonRecHitBuilder");

  theRPCInTheFit = par.getParameter<bool>("RefitRPCHits");

  theDYTthrs     = par.getParameter< std::vector<int> >("DYTthrs");
  theDYTselector = par.existsAs<int>("DYTselector")?par.getParameter<int>("DYTselector"):1;
  theDYTupdator = par.existsAs<bool>("DYTupdator")?par.getParameter<bool>("DYTupdator"):false;
  theDYTuseAPE = par.existsAs<bool>("DYTuseAPE")?par.getParameter<bool>("DYTuseAPE"):false;
  dytInfo        = new reco::DYTInfo();

  printStuff = par.getParameter<bool>("printStuff");


  if (par.existsAs<double>("RescaleErrorFactor")) {
    theRescaleErrorFactor = par.getParameter<double>("RescaleErrorFactor");
    edm::LogWarning("HighPTMuonRefitter") << "using error rescale factor " << theRescaleErrorFactor;
  }
  else theRescaleErrorFactor = 1000.;

  minNumHits = par.getParameter<int>("minNumHits");

  theCacheId_TRH = 0;
  theDTRecHitToken=iC.consumes<DTRecHitCollection>(theDTRecHitLabel);
  theCSCRecHitToken=iC.consumes<CSCRecHit2DCollection>(theCSCRecHitLabel);
  theGEMRecHitToken=iC.consumes<GEMRecHitCollection>(theGEMRecHitLabel);
  theME0RecHitToken=iC.consumes<ME0SegmentCollection>(theME0RecHitLabel); 
  CSCSegmentsToken = iC.consumes<CSCSegmentCollection>(InputTag("cscSegments"));
  all4DSegmentsToken=iC.consumes<DTRecSegment4DCollection>(InputTag("dt4DSegments"));

}

//--------------
// Destructor --
//--------------

HighPTMuonRefitter::~HighPTMuonRefitter() {
  delete dytInfo;
}


//
// set Event
//
void HighPTMuonRefitter::setEvent(const edm::Event& event) {

  theEvent = &event;
  event.getByToken(theDTRecHitToken, theDTRecHits);
  event.getByToken(theCSCRecHitToken, theCSCRecHits);
  event.getByToken(theGEMRecHitToken, theGEMRecHits);   
  event.getByToken(theME0RecHitToken, theME0RecHits);   
  event.getByToken(CSCSegmentsToken, CSCSegments);
  event.getByToken(all4DSegmentsToken, all4DSegments);
}


void HighPTMuonRefitter::setServices(const EventSetup& setup) {

  edm::ESHandle<TrajectoryFitter> aFitter;
  edm::ESHandle<TrajectorySmoother> aSmoother;
  setup.get<TrajectoryFitter::Record>().get(theFitterName,aFitter);
  setup.get<TrajectoryFitter::Record>().get(theSmootherName,aSmoother);
  theFitter = aFitter->clone();
  theSmoother.reset(aSmoother->clone());

  // Transient Rechit Builders
  unsigned long long newCacheId_TRH = setup.get<TransientRecHitRecord>().cacheIdentifier();
  if ( newCacheId_TRH != theCacheId_TRH ) {
    LogDebug(theCategory) << "TransientRecHitRecord changed!";
    setup.get<TransientRecHitRecord>().get(theTrackerRecHitBuilderName,theTrackerRecHitBuilder);
    setup.get<TransientRecHitRecord>().get(theMuonRecHitBuilderName,theMuonRecHitBuilder);
    hitCloner = static_cast<TkTransientTrackingRecHitBuilder const *>(theTrackerRecHitBuilder.product())->cloner();
  }
  theFitter->setHitCloner(&hitCloner);
  theSmoother->setHitCloner(&hitCloner);

}


//
// build a combined tracker-muon trajectory
//
pair<Trajectory,Trajectory>
HighPTMuonRefitter::refit(
		const reco::Track& globalTrack, 
		const std::string& theMuonHitsOption,
		const TrackerTopology *tTopo,
		const int& nHits,
		const int& hitPattern) const {
  LogTrace(theCategory) << " *** HighPTMuonRefitter *** option " << theMuonHitsOption << endl;
	/*
	std::cout << "\n\n************************************************" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << "************************************************" << std::endl;
	std::cout << " *** HighPTMuonRefitter *** option " << theMuonHitsOption << std::endl;
	std::cout << " *** HighPTMuonRefitter *** hitPattern " << hitPattern << std::endl;
	*/
    
	map<DetId, int> hitMap;
  // Get the TransientTrack
  reco::TransientTrack track(globalTrack,&*(theService->magneticField()),theService->trackingGeometry());
  
  // Get RecHits
  ConstRecHitContainer allRecHitsTemp = getRecHits(track);

  // Do tracker only
  if (theMuonHitsOption=="tracker") {
    ConstRecHitContainer trackerRecHits = removeMuonHits(allRecHitsTemp);
    vector<Trajectory> trackerOnlyTraj = transform(globalTrack,track,trackerRecHits);
    if (trackerOnlyTraj.empty()) {
      LogDebug(theCategory) << "No Tracker Track refitted!" << endl;
      return std::make_pair(Trajectory(),Trajectory());
    }
	
		vector<Trajectory> trackerTrajSM = theSmoother->trajectories(trackerOnlyTraj.front());
		if (!trackerTrajSM.empty()) return std::make_pair(trackerTrajSM.front(),Trajectory());
		else return std::make_pair(trackerOnlyTraj.front(),Trajectory());
  }

  // Select the Muon RecHits
  ConstRecHitContainer muonRecHitsTmp = applyMuonHitsOption(globalTrack,track,allRecHitsTemp,theMuonHitsOption,nHits,hitPattern,hitMap);
  // Remove Tracker RecHits
  ConstRecHitContainer muonRecHitsToRefit = removeTrackerHits(muonRecHitsTmp);

  // Obtain the Trajectory
  vector<Trajectory> muonOnlyTraj = transform(globalTrack,track,muonRecHitsToRefit);

  if (!muonOnlyTraj.size()) {
    LogTrace(theCategory) << "No refitted Tracks... " << endl;
    return std::make_pair(Trajectory(),Trajectory());
  } else {
    LogTrace(theCategory) << "Refitted pt: " 
			<< muonOnlyTraj.front().firstMeasurement().updatedState().globalParameters().momentum().perp() << endl;
  }

	// Update mu-only trajectory with vtx
	FreeTrajectoryState updatedFTSatVtx = 
		theHighPTUtilities->updateMuonOnlyTraj(muonOnlyTraj.front(), globalTrack);
	if (!updatedFTSatVtx.hasError()) return std::make_pair(Trajectory(),Trajectory());
	// Make trajectory
	vector<Trajectory> muonOnlyVtxTraj = 
		trajFromFTS(muonOnlyTraj.front(), updatedFTSatVtx, muonRecHitsToRefit);

	if (!muonOnlyVtxTraj.size()) {
		return std::make_pair(Trajectory(), Trajectory());
	}
	else {
		if (!muonOnlyVtxTraj.front().isValid() or !muonOnlyTraj.front().isValid()) {
			return std::make_pair(Trajectory(), Trajectory());
		}
		else {
			//std::cout << "After recalculating chi2" << std::endl;
			//std::cout << muonOnlyVtxTraj.front().chiSquared() << std::endl;
			//std::cout << std::endl;
			//printTrajectories(muonOnlyVtxTraj.front(),false);
			return std::make_pair(muonOnlyTraj.front(), muonOnlyVtxTraj.front());
		}
	}

}

//
// Get and return RecHits from TransientTrack
//
HighPTMuonRefitter::ConstRecHitContainer 
HighPTMuonRefitter::getRecHits(const reco::TransientTrack& track) const {
	auto tkbuilder = static_cast<TkTransientTrackingRecHitBuilder const *>(theTrackerRecHitBuilder.product());

	ConstRecHitContainer allRecHitsTemp;
	for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
		if ((*hit)->isValid()) {
			if ((*hit)->geographicalId().det() == DetId::Tracker) {
				allRecHitsTemp.push_back((**hit).cloneForFit(*tkbuilder->geometry()->idToDet( (**hit).geographicalId() ) ) );
			} // if hit is tracker
			else if ((*hit)->geographicalId().det() == DetId::Muon) {
				if ((*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit) {
					LogTrace(theCategory) << "RPC Rec Hit discarded"; 
					continue;
				} // if hit is RPC
				allRecHitsTemp.push_back(theMuonRecHitBuilder->build(&**hit));
			} // if hit is muon
		} // if hit is valid
	} // for loop over hits
	return allRecHitsTemp;
}

//
// Apply various RecHit selections
// default = all muon hits
//
HighPTMuonRefitter::ConstRecHitContainer 
HighPTMuonRefitter::applyMuonHitsOption(
		const reco::Track& globalTrack, const reco::TransientTrack& track,
		const ConstRecHitContainer& allRecHits, 
		const string& muonHitsOption, const int& nHits, const int& hitPattern,
		map<DetId, int> &hitMap) const {

	ConstRecHitContainer muonRecHitsForRefit;
	if (muonHitsOption=="default") {
		return allRecHits;
	}
	else if (muonHitsOption=="picky") {
		vector <Trajectory> globalTraj = transform(globalTrack, track, allRecHits);
		ConstRecHitContainer allRecHitsTmp = allRecHits;
		checkMuonHits(globalTrack, allRecHitsTmp, hitMap);
		muonRecHitsForRefit = selectMuonHits(globalTraj.front(),hitMap);
		return muonRecHitsForRefit;
	}
	else if (muonHitsOption=="firstHit") {
		muonRecHitsForRefit = getFirstHits(allRecHits);
		return muonRecHitsForRefit;
	}
	else if (muonHitsOption=="dyt") {
		vector <Trajectory> globalTraj = transform(globalTrack, track, allRecHits);
		//
		// DYT 2.0 
		//
		ConstRecHitContainer DYTRecHits;
		DynamicTruncation dytRefit(*theEvent,*theService);
		dytRefit.setProd(all4DSegments, CSCSegments);
		dytRefit.setSelector(theDYTselector);
		dytRefit.setThr(theDYTthrs);
		dytRefit.setUpdateState(theDYTupdator);
		dytRefit.setUseAPE(theDYTuseAPE);
		DYTRecHits = dytRefit.filter(globalTraj.front());
		dytInfo->CopyFrom(dytRefit.getDYTInfo());
		if ((DYTRecHits.size() > 1) && 
				(DYTRecHits.front()->globalPosition().mag() > DYTRecHits.back()->globalPosition().mag()))
			stable_sort(DYTRecHits.begin(),DYTRecHits.end(),RecHitLessByDet(alongMomentum));
		return DYTRecHits;
	}
	else if (muonHitsOption=="combinatoric") {
		ConstRecHitContainer muonRecHits = removeTrackerHits(allRecHits);
		muonRecHitsForRefit = applyHitMask(muonRecHits,nHits,hitPattern);
		return muonRecHitsForRefit;
		//
	}
	else {
		// not a valid muonHitsOption
		return muonRecHitsForRefit;
	}
}

//
// Apply combinatoric hit mask to muon rec hits
//
HighPTMuonRefitter::ConstRecHitContainer
HighPTMuonRefitter::applyHitMask(
		const ConstRecHitContainer &muonRecHits, const int& nHits, const int &hitMaskInt) const {
	
	ConstRecHitContainer selectedMuonRecHits;
	// Convert hit mask integer to binary
	std::vector<int> hitMask = convertToBinary(nHits, hitMaskInt);
	//std::cout << "*** before applying hit mask ***" << std::endl;
	//for (auto rh : muonRecHits) printRecHit(rh);
	//std::cout << "Hit mask" << std::endl;
	//for (auto h : hitMask) std::cout << h << " " ;
	//std::cout << std::endl;
	int ihit = 0;
	for (auto hit : muonRecHits) {
		// Least significant bit is the outermost Rec Hit in list
		if (hitMask[ihit]) {
			selectedMuonRecHits.push_back(hit);
		}
		ihit++;
	}
	/*
	std::cout << "*** after applying hit mask ***" << std::endl;
	for (auto rh : selectedMuonRecHits) printRecHit(rh);
	*/
	return selectedMuonRecHits;
}

//
// Convert hit mask integer (1 to 2**nhits-1) to binary
//
vector<int> HighPTMuonRefitter::convertToBinary(const int& nHits, const int &hitMaskInt) const {
	// Least significant bit is the outermost Rec Hit in list
	vector<int> hitMask(nHits,0);
	int num(hitMaskInt);
	int bin;
	int i=nHits-1;
	while (num>0) {
		bin = num % 2;
		hitMask[i] = bin;
		num /= 2;
		i--;
	}
	return hitMask;
}

//
// Remove Tracker Rec Hits
//
HighPTMuonRefitter::ConstRecHitContainer HighPTMuonRefitter::removeTrackerHits(const ConstRecHitContainer& hits) const
{
  ConstRecHitContainer results;
  ConstRecHitContainer::const_iterator it = hits.begin();
  for (; it!=hits.end(); it++) {

    DetId id = (*it)->geographicalId();

    //Check that this is a Muon hit that we're toying with -- else pass on this because the hacker is a moron / not careful

    if (id.det() == DetId::Tracker) {
		continue;
    }
    results.push_back(*it);
  }
  return results;
}

//
// Remove Muon Rec Hits
//
HighPTMuonRefitter::ConstRecHitContainer HighPTMuonRefitter::removeMuonHits(const ConstRecHitContainer& hits) const
{
  ConstRecHitContainer results;
  ConstRecHitContainer::const_iterator it = hits.begin();
  for (; it!=hits.end(); it++) {

    DetId id = (*it)->geographicalId();

    //Check that this is a Muon hit that we're toying with -- else pass on this because the hacker is a moron / not careful

    if (id.det() == DetId::Muon) {
			continue;
    }
    results.push_back(*it);
  }
  return results;
}

//
// Part of Picky
//
void HighPTMuonRefitter::checkMuonHits(const reco::Track& muon, 
				       ConstRecHitContainer& all,
				       map<DetId, int> &hitMap) const {

  LogTrace(theCategory) << " HighPTMuonRefitter::checkMuonHits " << endl;

  float coneSize = 20.0;

  // loop through all muon hits and calculate the maximum # of hits in each chamber
  for (ConstRecHitContainer::const_iterator imrh = all.begin(); imrh != all.end(); imrh++ ) {
        
    if ( (*imrh != 0 ) && !(*imrh)->isValid() ) continue;
  
    int detRecHits = 0;
    MuonRecHitContainer dRecHits;
      
    DetId id = (*imrh)->geographicalId();
    DetId chamberId;

    // Skip tracker hits
    if (id.det()!=DetId::Muon) continue;

    if ( id.subdetId() == MuonSubdetId::DT ) {
      DTChamberId did(id.rawId());
      chamberId=did;
      
      if ((*imrh)->recHits().size()>1) {
        std::vector <const TrackingRecHit*> hits2d = (*imrh)->recHits();
        for (std::vector <const TrackingRecHit*>::const_iterator hit2d = hits2d.begin(); hit2d!= hits2d.end(); hit2d++) {
          if ((*imrh)->recHits().size()>1) {
            std::vector <const TrackingRecHit*> hits1d = (*hit2d)->recHits();
            for (std::vector <const TrackingRecHit*>::const_iterator hit1d = hits1d.begin(); hit1d!= hits1d.end(); hit1d++) {
              DetId id1 = (*hit1d)->geographicalId();
              DTLayerId lid(id1.rawId());
              // Get the 1d DT RechHits from this layer
              DTRecHitCollection::range dRecHits = theDTRecHits->get(lid);
              int layerHits=0;
              for (DTRecHitCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
          	double rhitDistance = fabs(ir->localPosition().x()-(**hit1d).localPosition().x()); 
        	if ( rhitDistance < coneSize ) layerHits++;
                LogTrace(theCategory) << "       " << (ir)->localPosition() << "  " << (**hit1d).localPosition()
                     << " Distance: " << rhitDistance << " recHits: " << layerHits << "  SL: " << lid.superLayer() << endl;
              }
              if (layerHits>detRecHits) detRecHits=layerHits;
            }
          }
        }
      
      } else {
        DTLayerId lid(id.rawId());
    
        // Get the 1d DT RechHits from this layer
        DTRecHitCollection::range dRecHits = theDTRecHits->get(lid);

        for (DTRecHitCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
  	  double rhitDistance = fabs(ir->localPosition().x()-(**imrh).localPosition().x());
  	  if ( rhitDistance < coneSize ) detRecHits++;
          LogTrace(theCategory)	<< "       " << (ir)->localPosition() << "  " << (**imrh).localPosition()
               << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
        }
      }
    }// end of if DT
    else if ( id.subdetId() == MuonSubdetId::CSC ) {
      CSCDetId did(id.rawId());
      chamberId=did.chamberId();

      if ((*imrh)->recHits().size()>1) {
        std::vector <const TrackingRecHit*> hits2d = (*imrh)->recHits();
        for (std::vector <const TrackingRecHit*>::const_iterator hit2d = hits2d.begin(); hit2d!= hits2d.end(); hit2d++) {
          DetId id1 = (*hit2d)->geographicalId();
          CSCDetId lid(id1.rawId());
          
          // Get the CSC Rechits from this layer
          CSCRecHit2DCollection::range dRecHits = theCSCRecHits->get(lid);      
          int layerHits=0;

          for (CSCRecHit2DCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
    	    double rhitDistance = (ir->localPosition()-(**hit2d).localPosition()).mag();
  	    if ( rhitDistance < coneSize ) layerHits++;
            LogTrace(theCategory) << ir->localPosition() << "  " << (**hit2d).localPosition()
  	           << " Distance: " << rhitDistance << " recHits: " << layerHits << endl;
          }
          if (layerHits>detRecHits) detRecHits=layerHits;
        }
      } else {
        // Get the CSC Rechits from this layer
        CSCRecHit2DCollection::range dRecHits = theCSCRecHits->get(did);      

        for (CSCRecHit2DCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
  	  double rhitDistance = (ir->localPosition()-(**imrh).localPosition()).mag();
	  if ( rhitDistance < coneSize ) detRecHits++;
          LogTrace(theCategory) << ir->localPosition() << "  " << (**imrh).localPosition()
	         << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
        }
      }
    } //end of CSC if
    else if ( id.subdetId() == MuonSubdetId::GEM ) {
      GEMDetId did(id.rawId());
      chamberId=did.chamberId();

      if ((*imrh)->recHits().size()>1) {
        std::vector <const TrackingRecHit*> hits2d = (*imrh)->recHits();
        for (std::vector <const TrackingRecHit*>::const_iterator hit2d = hits2d.begin(); hit2d!= hits2d.end(); hit2d++) {
          DetId id1 = (*hit2d)->geographicalId();
          GEMDetId lid(id1.rawId());

          // Get the GEM Rechits from this layer
          GEMRecHitCollection::range dRecHits = theGEMRecHits->get(lid);
          int layerHits=0;

          for (GEMRecHitCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
            double rhitDistance = (ir->localPosition()-(**hit2d).localPosition()).mag();
            if ( rhitDistance < coneSize ) layerHits++;
            LogTrace(theCategory) << ir->localPosition() << "  " << (**hit2d).localPosition()
                   << " Distance: " << rhitDistance << " recHits: " << layerHits << endl;
          }
          if (layerHits>detRecHits) detRecHits=layerHits;
        }
      } else {
        // Get the GEM Rechits from this layer
        GEMRecHitCollection::range dRecHits = theGEMRecHits->get(did);

        for (GEMRecHitCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
          double rhitDistance = (ir->localPosition()-(**imrh).localPosition()).mag();
          if ( rhitDistance < coneSize ) detRecHits++;
          LogTrace(theCategory) << ir->localPosition() << "  " << (**imrh).localPosition()
                 << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
        }
      }
    } //end of GEM if
    else if ( id.subdetId() == MuonSubdetId::ME0 ) {
      ME0DetId did(id.rawId());
      chamberId=did.chamberId();

      if ((*imrh)->recHits().size()>1) {
        std::vector <const TrackingRecHit*> hits2d = (*imrh)->recHits();
        for (std::vector <const TrackingRecHit*>::const_iterator hit2d = hits2d.begin(); hit2d!= hits2d.end(); hit2d++) {
          DetId id1 = (*hit2d)->geographicalId();
          ME0DetId lid(id1.rawId());

          // Get the ME0 Rechits from this layer
          ME0SegmentCollection::range dRecHits = theME0RecHits->get(lid);
          int layerHits=0;

          for (ME0SegmentCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
            double rhitDistance = (ir->localPosition()-(**hit2d).localPosition()).mag();
            if ( rhitDistance < coneSize ) layerHits++;
            LogTrace(theCategory) << ir->localPosition() << "  " << (**hit2d).localPosition()
                   << " Distance: " << rhitDistance << " recHits: " << layerHits << endl;
          }
          if (layerHits>detRecHits) detRecHits=layerHits;
        }
      } else {
        // Get the ME0 Rechits from this layer
        ME0SegmentCollection::range dRecHits = theME0RecHits->get(did);

        for (ME0SegmentCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
          double rhitDistance = (ir->localPosition()-(**imrh).localPosition()).mag();
          if ( rhitDistance < coneSize ) detRecHits++;
          LogTrace(theCategory) << ir->localPosition() << "  " << (**imrh).localPosition()
                 << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
        }
      }
    } //end of ME0 if
    else {
      if ( id.subdetId() != MuonSubdetId::RPC ) LogError(theCategory)<<" Wrong Hit Type ";
      continue;      
    }
      
    map<DetId,int>::iterator imap=hitMap.find(chamberId);
    if (imap!=hitMap.end()) {
      if (detRecHits>imap->second) imap->second=detRecHits;
    } else hitMap[chamberId]=detRecHits;

  } // end of loop over muon rechits

  for (map<DetId,int>::iterator imap=hitMap.begin(); imap!=hitMap.end(); imap++ ) 
    LogTrace(theCategory) << " Station " << imap->first.rawId() << ": " << imap->second <<endl; 

  LogTrace(theCategory) << "CheckMuonHits: "<<all.size();

  // check order of muon measurements
  if ( (all.size() > 1) &&
       ( all.front()->globalPosition().mag() >
	 all.back()->globalPosition().mag() ) ) {
    LogTrace(theCategory)<< "reverse order: ";
    stable_sort(all.begin(),all.end(),RecHitLessByDet(alongMomentum));
  }
}


//
// Get the hits from the first muon station (containing hits)
//
HighPTMuonRefitter::ConstRecHitContainer 
HighPTMuonRefitter::getFirstHits(const ConstRecHitContainer& all) const {

  LogTrace(theCategory) << " HighPTMuonRefitter::getFirstHits\nall rechits length:" << all.size() << endl;
  ConstRecHitContainer first;

  int station_to_keep = 999;
  vector<int> stations;
  for (ConstRecHitContainer::const_iterator ihit = all.begin(); ihit != all.end(); ++ihit) {
  
    int station = 0;
    bool use_it = true;
    DetId id = (*ihit)->geographicalId();
    unsigned raw_id = id.rawId();
    if (!(*ihit)->isValid()) station = -1;
          else {
	if (id.det() == DetId::Muon) {
	  switch (id.subdetId()) {
	  case MuonSubdetId::DT:  station = DTChamberId(raw_id).station(); break;
	  case MuonSubdetId::CSC: station = CSCDetId(raw_id).station(); break;
	  case MuonSubdetId::GEM: station = GEMDetId(raw_id).station(); break;
	  case MuonSubdetId::ME0: station = ME0DetId(raw_id).station(); break;
	  case MuonSubdetId::RPC: station = RPCDetId(raw_id).station(); use_it = false; break;
	  }
	}
      }


    if (use_it && station > 0 && station < station_to_keep) station_to_keep = station;
    stations.push_back(station);
    LogTrace(theCategory) << "rawId: " << raw_id << " station = " << station << " station_to_keep is now " << station_to_keep;
  }

  if (station_to_keep <= 0 || station_to_keep > 4 || stations.size() != all.size())
    LogInfo(theCategory) << "failed to getFirstHits (all muon hits are outliers/bad ?)! station_to_keep = " 
			    << station_to_keep << " stations.size " << stations.size() << " all.size " << all.size();

  for (unsigned i = 0; i < stations.size(); ++i)
    if (stations[i] >= 0 && stations[i] <= station_to_keep) first.push_back(all[i]);

  return first;
}


//
// Picky
// select muon hits compatible with trajectory; 
// check hits in chambers with showers
//
HighPTMuonRefitter::ConstRecHitContainer 
HighPTMuonRefitter::selectMuonHits(const Trajectory& traj, 
                                   const map<DetId, int> &hitMap) const {

  ConstRecHitContainer muonRecHits;
  const double globalChi2Cut = 200.0;

  vector<TrajectoryMeasurement> muonMeasurements = traj.measurements(); 

  // loop through all muon hits and skip hits with bad chi2 in chambers with high occupancy      
  for (std::vector<TrajectoryMeasurement>::const_iterator im = muonMeasurements.begin(); im != muonMeasurements.end(); im++ ) {

    if ( !(*im).recHit()->isValid() ) continue;
    if ( (*im).recHit()->det()->geographicalId().det() != DetId::Muon ) {
      //      if ( ( chi2ndf < globalChi2Cut ) )
      muonRecHits.push_back((*im).recHit());
      continue;
    }  
    const MuonTransientTrackingRecHit* immrh = dynamic_cast<const MuonTransientTrackingRecHit*>((*im).recHit().get());

    DetId id = immrh->geographicalId();
    DetId chamberId;
    int threshold = 0;
    double chi2Cut = 0.0;

    // get station of hit if it is in DT
    if ( (*immrh).isDT() ) {
      DTChamberId did(id.rawId());
      chamberId = did;
      threshold = theHitThreshold;
      chi2Cut = theDTChi2Cut;
    }
    // get station of hit if it is in CSC
    else if ( (*immrh).isCSC() ) {
      CSCDetId did(id.rawId());
      chamberId = did.chamberId();
      threshold = theHitThreshold;
      chi2Cut = theCSCChi2Cut;
    }
    // get station of hit if it is in GEM
    else if ( (*immrh).isGEM() ) {
      GEMDetId did(id.rawId());
      chamberId = did.chamberId();
      threshold = theHitThreshold;
      chi2Cut = theGEMChi2Cut;
    }
    // get station of hit if it is in ME0
    else if ( (*immrh).isME0() ) {
      ME0DetId did(id.rawId());
      chamberId = did.chamberId();
      threshold = theHitThreshold;
      chi2Cut = theME0Chi2Cut;
    }
    // get station of hit if it is in RPC
    else if ( (*immrh).isRPC() ) {
      RPCDetId rpcid(id.rawId());
      chamberId = rpcid;
      threshold = theHitThreshold;
      chi2Cut = theRPCChi2Cut;
    } else
      continue;

    double chi2ndf = (*im).estimate()/(*im).recHit()->dimension();  

    bool keep = true;
    map<DetId,int>::const_iterator imap=hitMap.find(chamberId);
    if ( imap!=hitMap.end() ) 
      if (imap->second>threshold) keep = false;
    
    if ( (keep || (chi2ndf<chi2Cut)) && (chi2ndf<globalChi2Cut) ) {
      muonRecHits.push_back((*im).recHit());
    } else {
      LogTrace(theCategory)
	<< "Skip hit: " << id.rawId() << " chi2=" 
	<< chi2ndf << " ( threshold: " << chi2Cut << ") Det: " 
	<< imap->second << endl;
    }
  }
  
  // check order of rechits
  reverse(muonRecHits.begin(),muonRecHits.end());
  return muonRecHits;
}




/*
double HighPTMuonRefitter::computeChi2(const Trajectory& traj) const {
	
	// Remove the RecHit that has the worst chi^2 wrt the N-1 track fit

  //std::cout << "\n================================================\n" << std::endl;
	//std::cout << "RecHit chi2 / ndf for muon hits\n" << std::endl;

	
  TrajectoryStateCombiner combiner;
  MeasurementEstimator *theEstimator;
  theEstimator = new Chi2MeasurementEstimator(100.,10.);
  ConstRecHitContainer muonRecHits;
  //const double globalChi2Cut = 200.0;

	double fullChi2(0.);
  vector<TrajectoryMeasurement> muonMeasurements = traj.measurements(); 

  int hitCounter = muonMeasurements.size();
  //int theWorstHit = -1;
  //double theWorstEstimate = -1;

  // loop through all muon hits and skip hits with bad chi2 in chambers with high occupancy      
  for (std::vector<TrajectoryMeasurement>::const_iterator im = muonMeasurements.begin(); im != muonMeasurements.end(); im++,--hitCounter ) {

	  // only include valid hits
    if ( !(*im).recHit()->isValid() ) continue;
	  // Add all tracker hits
    if ( (*im).recHit()->det()->geographicalId().det() != DetId::Muon ) {
      //      if ( ( chi2ndf < globalChi2Cut ) )
      muonRecHits.push_back((*im).recHit());
      continue;
    }  

    const MuonTransientTrackingRecHit* immrh = dynamic_cast<const MuonTransientTrackingRecHit*>((*im).recHit().get());

	TrajectoryStateOnSurface combFwdBwd;
	// For last hit, smoothed prediction is forward prediction
	// For all other hits, smoothed prediction is combination of forward and backward predictions
	if (hitCounter==1) {
		combFwdBwd = (*im).forwardPredictedState();
	} else {
		TrajectoryStateOnSurface fwdPred = (*im).forwardPredictedState();
		TrajectoryStateOnSurface bwdPred = (*im).backwardPredictedState();
		combFwdBwd = combiner(fwdPred,bwdPred);
	}
	
	double thisEstimate;
	thisEstimate = theEstimator->estimate(combFwdBwd, *immrh).second;
	fullChi2 += thisEstimate;


	//std::cout << "Combined Fwd+Bwd TSOS" << std::endl;
	//printTSOS(combFwdBwd);
  //std::cout << "\n------------------------------------------------\n" << std::endl;
	//printRecHit((*im).recHitP());
  //std::cout << "\n------------------------------------------------\n" << std::endl;
	//std::cout << (*immrh).det()->subDetector() << " " << thisEstimate << " " << (*immrh).dimension() << " " << thisEstimate/(*immrh).dimension() << std::endl;
  //std::cout << "\n================================================\n" << std::endl;

	if (thisEstimate>theWorstEstimate) {
		theWorstEstimate = thisEstimate;
		theWorstHit = thisHit;
	}
 // muonRecHits.push_back((*im).recHit());



  }
	//std::cout << std::endl;
  //return muonRecHits;
	return fullChi2;
}
*/





//
// Determine RecHit ordering (inside-out, etc.)
//
HighPTMuonRefitter::RefitDirection
HighPTMuonRefitter::checkRecHitsOrdering(const TransientTrackingRecHit::ConstRecHitContainer& recHits) const {

  if (!recHits.empty()){
    ConstRecHitContainer::const_iterator frontHit = recHits.begin();
    ConstRecHitContainer::const_iterator backHit  = recHits.end() - 1;
    while( !(*frontHit)->isValid() && frontHit != backHit) {frontHit++;}
    while( !(*backHit)->isValid() && backHit != frontHit)  {backHit--;}

    double rFirst = (*frontHit)->globalPosition().mag();
    double rLast  = (*backHit) ->globalPosition().mag();

    if(rFirst < rLast) return insideOut;
    else if(rFirst > rLast) return outsideIn;
    else {
      LogError(theCategory) << "Impossible determine the rechits order" <<endl;
      return undetermined;
    }
  } else {
    LogError(theCategory) << "Impossible determine the rechits order" <<endl;
    return undetermined;
  }
}


//
// Convert Tracks into Trajectories with a given set of hits
//
vector<Trajectory> HighPTMuonRefitter::transform(const reco::Track& newTrack,
						 const reco::TransientTrack track,
						 const TransientTrackingRecHit::ConstRecHitContainer& urecHitsForReFit) const {
  
	TransientTrackingRecHit::ConstRecHitContainer recHitsForReFit = urecHitsForReFit;
	LogTrace(theCategory) << "HighPTMuonRefitter::transform: " << recHitsForReFit.size() << " hits:";

	if(recHitsForReFit.size() < minNumHits) return vector<Trajectory>();

	// Check the order of the rechits
	RefitDirection recHitsOrder;
	if (recHitsForReFit.size()==1) recHitsOrder = theRefitDirection;
	else recHitsOrder = checkRecHitsOrdering(recHitsForReFit);

	LogTrace(theCategory) << "checkRecHitsOrdering() returned " << recHitsOrder
		<< ", theRefitDirection is " << theRefitDirection
		<< " (insideOut == " << insideOut << ", outsideIn == " << outsideIn << ")";

	// Reverse the order in the case of inconsistency between the fit direction and the rechit order
	if(theRefitDirection != recHitsOrder) reverse(recHitsForReFit.begin(),recHitsForReFit.end());

	// Even though we checked the rechits' ordering above, we may have
	// already flipped them elsewhere (getFirstHits() is such a
	// culprit). Use the global positions of the states and the desired
	// refit direction to find the starting TSOS.
	TrajectoryStateOnSurface firstTSOS, lastTSOS;
	unsigned int innerId; //UNUSED: outerId;
	bool order_swapped = track.outermostMeasurementState().globalPosition().mag() < track.innermostMeasurementState().globalPosition().mag();
	bool inner_is_first;
	LogTrace(theCategory) << "order swapped? " << order_swapped;

	// Fill the starting state, depending on the ordering above.
	if (printStuff) {
		if (order_swapped) std::cout << " order swapped : " << order_swapped << std::endl;
		std::cout << " checkRecHitsOrdering() returned " << recHitsOrder << std::endl;
		std::cout << " theRefitDirection is " << theRefitDirection << std::endl;
		std::cout << " (insideOut==" << insideOut << ", outsideIn==" << outsideIn << ")" << std::endl;
	}
	if ( (theRefitDirection == insideOut && !order_swapped) 
			|| (theRefitDirection == outsideIn && order_swapped)) {
		if (printStuff) {
			std::cout << " innerId = newTrack.innerDetId()" << std::endl;
			std::cout << " firstTSOS = track.innermostMeasurementstate()" <<std::endl;
			std::cout << " lastTSOS = track.outermostMeasurementstate()" <<std::endl;
			std::cout << " inner_is_first = true" <<std::endl;
		}
		innerId   = newTrack.innerDetId();
		//UNUSED:    outerId   = newTrack.outerDetId();
		firstTSOS = track.innermostMeasurementState();
		lastTSOS  = track.outermostMeasurementState();
		inner_is_first = true;
	}
	else {
		if (printStuff) {
			std::cout << " innerId = newTrack.outerDetId()" << std::endl;
			std::cout << " firstTSOS = track.outermostMeasurementstate()" <<std::endl;
			std::cout << " lastTSOS = track.innermostMeasurementstate()" <<std::endl;
			std::cout << " inner_is_first = false" <<std::endl;
		}
		innerId   = newTrack.outerDetId();
		//UNUSED:    outerId   = newTrack.innerDetId();
		firstTSOS = track.outermostMeasurementState();
		lastTSOS  = track.innermostMeasurementState();
		inner_is_first = false;
	} 

	if (printStuff) {
		if (inner_is_first) std::cout << "firstTSOS: inner_is_first : " << inner_is_first << std::endl;
		std::cout << " globalPosition is " << firstTSOS.globalPosition() << std::endl;
		std::cout << " innerId is " << innerId << std::endl;
	}

	LogTrace(theCategory) << "firstTSOS: inner_is_first? " << inner_is_first
		<< " globalPosition is " << firstTSOS.globalPosition()
		<< " innerId is " << innerId;

	if(!firstTSOS.isValid()){
		LogWarning(theCategory) << "Error wrong initial state!" << endl;
		if (printStuff) std::cout << "Error wrong initial state!" << std::endl;
		return vector<Trajectory>();
	}
  
	if (printStuff) {
		std::cout << "\n------------------------------------------------" << std::endl;
		std::cout << "\nFirst TSOS before rescaling error" << std::endl;
		printTSOS(firstTSOS);
	}


	// This is the only way to get a TrajectorySeed with settable propagation direction
	PTrajectoryStateOnDet garbage1;
	edm::OwnVector<TrackingRecHit> garbage2;
	PropagationDirection propDir;

	// These lines cause the code to ignore completely what was set
	// above, and force propDir for tracks from collisions!
	if(propDir == alongMomentum && theRefitDirection == outsideIn)  propDir=oppositeToMomentum;
	if(propDir == oppositeToMomentum && theRefitDirection == insideOut) propDir=alongMomentum;

	// For now, don't put any special propagation direction logic for cosmic muons
	// TrackingTools/TrackRefitters/src/TrackTransformer.cc has correct logic to handle this

	TrajectorySeed seed(garbage1,garbage2,propDir);

	if(recHitsForReFit.front()->geographicalId() != DetId(innerId)){
		LogDebug(theCategory)<<"Propagation occured"<<endl;
		if (printStuff) std::cout<<"Propagation occured"<<std::endl;

		LogTrace(theCategory) << "propagating firstTSOS at " << firstTSOS.globalPosition()
			  << " to first rechit with surface pos " 
			  << recHitsForReFit.front()->det()->surface().toGlobal(LocalPoint(0,0,0));
		if (printStuff) {
			std::cout << "propagating firstTSOS at " << firstTSOS.globalPosition()
					  << " to first rechit with surface pos " 
					  << recHitsForReFit.front()->det()->surface().toGlobal(LocalPoint(0,0,0)) 
					  << std::endl;
		}

		firstTSOS = theService->propagator(thePropagatorName)->propagate(firstTSOS, recHitsForReFit.front()->det()->surface());
		if(!firstTSOS.isValid()){
			LogDebug(theCategory)<<"Propagation error!"<<endl;
			if (printStuff) std::cout<<"Propagation error!"<<std::endl;
			return vector<Trajectory>();
		}
	}
	if (printStuff) {
		std::cout << "\n------------------------------------------------" << std::endl;
		std::cout << "First TSOS after propagation to surface of first RecHit" << std::endl;
		printTSOS(firstTSOS);
	}

	firstTSOS.rescaleError(theRescaleErrorFactor);

	if (printStuff) {
		std::cout << "\n------------------------------------------------" << std::endl;
		std::cout << "First TSOS after rescaling error" << std::endl;
		printTSOS(firstTSOS);
	}

       
	if(!firstTSOS.isValid()){ std::cout << "first TSOS is not valid?" << std::endl;}

	vector<Trajectory> trajectories = theFitter->fit(seed,recHitsForReFit,firstTSOS);

	if(trajectories.empty()){
		LogDebug(theCategory) << "No Track refitted!" << endl;
		if (printStuff) std::cout << "No Track refitted!" << std::endl;
		return vector<Trajectory>();
	}

	return trajectories;
}

//
// Build trajectory from FTS at vertex and rechits used in refit
// Calculate chi^2 propagating inside-out from the FTS at vertex to each RH 
//
vector<Trajectory>
HighPTMuonRefitter::trajFromFTS(
		const Trajectory& muonOnlyTraj,
		const FreeTrajectoryState& ftsAtVtx, const ConstRecHitContainer& recHits) const {

	// Need a non-const version in case order needs to change
	TransientTrackingRecHit::ConstRecHitContainer recHitsToUse = recHits;
	
	// Make trajectory
	PropagationDirection propDir = alongMomentum;
	PTrajectoryStateOnDet garbage1;
	edm::OwnVector<TrackingRecHit> garbage2;
	TrajectorySeed seed(garbage1,garbage2,propDir);
	Trajectory traj(seed,propDir);

	// Set Estimator
  MeasurementEstimator *theEstimator;
  theEstimator = new Chi2MeasurementEstimator(100.,10.);


	// Force the RH ordering to be inside-out
	RefitDirection recHitsOrder;
	if (recHitsToUse.size()>1) recHitsOrder = checkRecHitsOrdering(recHitsToUse);
	else recHitsOrder = insideOut;

	if (recHitsOrder==outsideIn){
		reverse(recHitsToUse.begin(), recHitsToUse.end());
	}

	// 
	// Propagate FTS to first RecHit surface
	//
	TrajectoryStateOnSurface firstTSOS = 
		theService->propagator(thePropagatorName)->propagate(ftsAtVtx, recHitsToUse.front()->det()->surface());
		//theService->propagator(thePropagatorName)->propagate(ftsAtVtx, *((*recHitsToUse.front()).surface()));
	if (!firstTSOS.isValid()) {
		return std::vector<Trajectory>();
	}
	TrajectoryStateOnSurface firstPredTSOS(firstTSOS);
	//
	//
	TrajectoryStateOnSurface thisTSOS, prevTSOS;
	int i(0);
  for (auto ihit : recHitsToUse) {

		const TransientTrackingRecHit & hit = (*ihit);

		if (i==0) {
			thisTSOS=firstPredTSOS;
		}
		else {
			//thisTSOS = theService->propagator(thePropagatorName)->propagate(prevTSOS, *(hit.surface()));
			thisTSOS = theService->propagator(thePropagatorName)->propagate(prevTSOS, ihit->det()->surface());
		}

		if (thisTSOS.isValid()){

			std::pair<bool,double> thisEstimate = theEstimator->estimate(thisTSOS, hit);
			if (!thisEstimate.first) {
				//std::cout << "HighPTMuonRefitter trajFromFTS thisEstimate.first is false" << std::endl;
				return std::vector<Trajectory>();
			}

			TrajectoryMeasurement thisMuonOnlyTM;
			for (auto tmpMuonOnlyTM : muonOnlyTraj.measurements()) {
				if (tmpMuonOnlyTM.recHit()->rawId()==hit.rawId()) {
					thisMuonOnlyTM=tmpMuonOnlyTM;
					break;
				}
			}

			TrajectoryMeasurement thisTM = 
				TrajectoryMeasurement(thisTSOS, ihit, thisEstimate.second, thisMuonOnlyTM.layer());

			traj.push(thisTM);
			prevTSOS = thisTSOS;
			i++;
		}
		else {
			//std::cout << "HighPTMuonRefitter trajFromFTS() thisTSOS is not valid" << std::endl;
			return std::vector<Trajectory>();
		}
	}
	return std::vector<Trajectory>(1,std::move(traj));
}


/*
 *
 * Printing Functions
 * (should eventually convert all this to LogInfo() or LogTrace() or whatever)
 *
 *
 */


//
// Print all info about trajectories
//
void HighPTMuonRefitter::printTrajectories(const Trajectory& trajectories, const bool &smoothing) const {
  std::cout << "Number of trajectory measurements " << trajectories.measurements().size() << std::endl;
	
  int hitcounter(0);
  if (smoothing) hitcounter = trajectories.measurements().size();
  for (auto tm : trajectories.measurements()) {
    if (!smoothing) hitcounter++;
	std::cout << "\nHit " << hitcounter << endl;

	std::cout << "\nForward Predicted State" << std::endl;
	printTSOS(tm.forwardPredictedState());
	if (smoothing) {
		std::cout << "\nBackward Predicted State" << std::endl;
		printTSOS(tm.backwardPredictedState());
	}

	std::cout << "\n------------------------------------------------" << std::endl;
	printRecHit(tm.recHitP());
	std::cout << "\n------------------------------------------------" << std::endl;
	std::cout << "\nRecHit chi2 w.r.t. Forward Predicted State" << std::endl;
    double chi2ndf = tm.estimate()/tm.recHit()->dimension();  
	std::cout << tm.estimate() << " / " << tm.recHit()->dimension() << " = " << chi2ndf << std::endl;
	std::cout << "\n------------------------------------------------" << std::endl;

	if (smoothing) {
		std::cout << "\nSmoothed State (fwd + bwd + hit)\n" << std::endl;
		printTSOS(tm.updatedState());
	} else {
		std::cout << "\nUpdated State (fwd + hit)\n" << std::endl;
		printTSOS(tm.updatedState());
	}
	if (smoothing) hitcounter--;
	std::cout << "\n================================================" << std::endl;
  }
  std::cout << "************************************************" << std::endl;
}
		
void HighPTMuonRefitter::printTSOS(const TrajectoryStateOnSurface& tsos) const {
	CurvilinearTrajectoryParameters curvil(tsos.globalPosition(),
			tsos.globalMomentum(), tsos.charge());
  std::cout
	  << "\nCurviliner Parameters\n" <<
	  curvil.vector()
	  << "\nCurvilinear Error\n" <<
	  tsos.curvilinearError().matrix()
	  << "\n\nGlobal position\n" << 
	  tsos.globalPosition()
	  << "\nGlobal Momemtum\n" << 
	  tsos.globalMomentum()
	  << "\n\npt = " << std::sqrt(tsos.globalMomentum().x()*tsos.globalMomentum().x() + tsos.globalMomentum().y()*tsos.globalMomentum().y()) 
	  << "\n\nGlobal Direction\n" <<
	  tsos.globalDirection()
	  << "\nCartesian Error\n" << 
	  tsos.cartesianError().matrix()
	  << "\n\nLocal Parameters\n" <<
	  tsos.localParameters().vector()
	  << "\nLocal Error\n" <<
	  tsos.localError().matrix()
	  << "\nLocal Position\n" <<
	  tsos.localPosition()
	  << "\nLocal Momentum\n" <<
	  tsos.localMomentum()
	  << "\nLocal Direction\n" <<
	  tsos.localDirection()
  << std::endl;
}

//
// print RecHits
//
void HighPTMuonRefitter::printHitscout(const ConstRecHitContainer& hits) const {

	std::cout << "\n================================================" << std::endl;
	std::cout << "\nList of All RecHits used in Track" << std::endl;
	std::cout << "\nUsed RecHits: " << hits.size() << std::endl;
  for (ConstRecHitContainer::const_iterator ir = hits.begin(); ir != hits.end(); ir++ ) {
	  auto hit = (*ir);
	  printRecHit(hit);
  std::cout << "\n------------------------------------------------\n" << std::endl;
  }
  std::cout << "\n================================================\n" << std::endl;

}

//
// print RecHits
//
void HighPTMuonRefitter::printHits(const ConstRecHitContainer& hits) const {

  LogTrace(theCategory) << "Used RecHits: " << hits.size();
  for (ConstRecHitContainer::const_iterator ir = hits.begin(); ir != hits.end(); ir++ ) {
    if ( !(*ir)->isValid() ) {
      LogTrace(theCategory) << "invalid RecHit";
      continue; 
    }
    
    const GlobalPoint& pos = (*ir)->globalPosition();
    
    LogTrace(theCategory) 
      << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
      << "  z = " << pos.z()
      << "  dimension = " << (*ir)->dimension()
      << "  det = " << (*ir)->det()->geographicalId().det()
      << "  subdet = " << (*ir)->det()->subDetector()
      << "  raw id = " << (*ir)->det()->geographicalId().rawId();
  }

}
void HighPTMuonRefitter::printRecHit(const TrackingRecHit::ConstRecHitPointer & hit) const {
	std::cout << "\nRecHit" << std::endl;
	if (hit->isValid()) std::cout << "\nHit is valid" << std::endl;
	else std::cout << "\nHit is invalid" << std::endl;
	const GlobalPoint& pos = hit->globalPosition();
	std::cout <<
	  "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y()) << "\n" <<
	  //"x, y, z = " << pos.x() << ", " << pos.y() << ", " << pos.z() << "\n" <<
	  "dimension = " << hit->dimension() << "\n" <<
	  "det = " << hit->det()->geographicalId().det() << "\n" <<
	  "subdet = " << hit->det()->subDetector() << "\n" <<
	  "raw id = " << hit->det()->geographicalId().rawId() <<
	std::endl;
	std::cout << "Projection Matrix" << std::endl;
	std::cout << hit->projectionMatrix() << std::endl;
	std::cout << "Parameters" << std::endl;
	std::cout << hit->parameters() << std::endl;
	std::cout << "Parameters Error" << std::endl;
	std::cout << hit->parametersError() << std::endl;
	std::cout << "Global Position" << std::endl;
	std::cout << hit->globalPosition() << std::endl;
	std::cout << "Global Position Error" << std::endl;
	std::cout << hit->globalPositionError().matrix() << std::endl;
	std::cout << "Local Position" << std::endl;
	std::cout << hit->localPosition() << std::endl;
	std::cout << "Local Position Error" << std::endl;
	std::cout << hit->localPositionError() << std::endl;
	/*
	if (hit->dimension()==4){
	std::cout << "Local Direction" << std::endl;
	std::cout << hit->localDirection() << std::endl;
	std::cout << "Local Direction Error" << std::endl;
	std::cout << hit->localDirectionError() << std::endl;
	std::cout << "chi2 / dof" << std::endl;
	std::cout << hit->chi2() << " / " << hit->degreesOfFreedom() << std::endl;
	}
	*/
}
