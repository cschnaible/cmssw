/**  \class HighPTMuonProducer
 * 
 *   TeV muon reconstructor:
 *
 *
 *
 *   \author Christian Schnaible (UCLA)
 */

// Framework
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "RecoMuon/GlobalMuonProducer/src/HighPTMuonProducer.h"

// TrackFinder and specific GLB Trajectory Builder
#include "RecoMuon/GlobalTrackFinder/interface/GlobalMuonTrajectoryBuilder.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackFinder.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackLoader.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonRefitter.h"
#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonUtilities.h"


using namespace edm;
using namespace std;
using namespace reco;

//
// constructor with config
//
HighPTMuonProducer::HighPTMuonProducer(const ParameterSet& parameterSet) {

  LogDebug("Muon|RecoMuon|HighPTMuonProducer") << "constructor called" << endl;

  // GLB Muon Collection Label
  theGLBCollectionLabel = parameterSet.getParameter<InputTag>("MuonCollectionLabel");
  glbMuonsToken=consumes<reco::TrackCollection>(theGLBCollectionLabel);
  glbMuonsTrajToken=consumes<std::vector<Trajectory> >(theGLBCollectionLabel.label());

  // service parameters
  ParameterSet serviceParameters = parameterSet.getParameter<ParameterSet>("ServiceParameters");

  // the services
  theService = new MuonServiceProxy(serviceParameters);
  edm::ConsumesCollector iC  = consumesCollector();  

	// Utilities parameters
  ParameterSet updatorParameters = parameterSet.getParameter<ParameterSet>("UtilitiesParameters");
	theHighPTUtilities = new HighPTMuonUtilities(updatorParameters, theService, iC);

  // TrackRefitter parameters
  ParameterSet refitterParameters = parameterSet.getParameter<ParameterSet>("RefitterParameters");
  theRefitter = new HighPTMuonRefitter(refitterParameters, theService, theHighPTUtilities, iC);
	
  // TrackLoader parameters
  ParameterSet trackLoaderParameters = parameterSet.getParameter<ParameterSet>("TrackLoaderParameters");
  theTrackLoader = new MuonTrackLoader(trackLoaderParameters,iC,theService);

  theRefits = parameterSet.getParameter< std::vector<std::string> >("Refits");

  for(unsigned int ww=0;ww<theRefits.size();ww++){
    LogDebug("Muon|RecoMuon|HighPTMuonProducer") << "Refit " << theRefits[ww];
    produces<reco::TrackCollection>(theRefits[ww]);
    produces<TrackingRecHitCollection>(theRefits[ww]);
    produces<reco::TrackExtraCollection>(theRefits[ww]);
    produces<vector<Trajectory> >(theRefits[ww]) ;
    produces<TrajTrackAssociationCollection>(theRefits[ww]);
	produces<reco::TrackToTrackMap>(theRefits[ww]);
	if (theRefits[ww]!="tracker"){
		LogDebug("Muon|RecoMuon|HighPTMuonProducer") << "Refit " << theRefits[ww]<<"VtxUpdate";
		produces<reco::TrackCollection>(theRefits[ww]+"VtxUpdate");
		produces<TrackingRecHitCollection>(theRefits[ww]+"VtxUpdate");
		produces<reco::TrackExtraCollection>(theRefits[ww]+"VtxUpdate");
		produces<vector<Trajectory> >(theRefits[ww]+"VtxUpdate") ;
		produces<TrajTrackAssociationCollection>(theRefits[ww]+"VtxUpdate");
		produces<reco::TrackToTrackMap>(theRefits[ww]+"VtxUpdate");
	}
  }
  produces<DYTestimators> ("dytInfo");
}


//
// destructor
//
HighPTMuonProducer::~HighPTMuonProducer() {

  LogTrace("Muon|RecoMuon|HighPTMuonProducer") << "destructor called" << endl;
  if (theService) delete theService;
  if (theRefitter) delete theRefitter;
  if (theTrackLoader) delete theTrackLoader;
	if (theHighPTUtilities) delete theHighPTUtilities;
}


//
// reconstruct muons
//
void HighPTMuonProducer::produce(Event& event, const EventSetup& eventSetup) {

  const string metname = "Muon|RecoMuon|HighPTMuonProducer";  
  LogTrace(metname)<< endl << endl;
  LogTrace(metname)<< "High pT Muon Reconstruction started" << endl;  

  // Update the services
  theService->update(eventSetup);

  theRefitter->setEvent(event);
  theHighPTUtilities->setEvent(event);

  theRefitter->setServices(theService->eventSetup());

  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoHand;
  eventSetup.get<TrackerTopologyRcd>().get(tTopoHand);
  const TrackerTopology *tTopo=tTopoHand.product();


  // Take the GLB muon container(s)
  Handle<reco::TrackCollection> glbMuons;
  event.getByToken(glbMuonsToken,glbMuons);

  auto dytInfo = std::make_unique<DYTestimators>();
  DYTestimators::Filler filler(*dytInfo);
  size_t GLBmuonSize = glbMuons->size();
  vector<DYTInfo> dytTmp(GLBmuonSize);

  Handle<vector<Trajectory> > glbMuonsTraj;

  LogTrace(metname)<< "Taking " << glbMuons->size() << " Global Muons "<<theGLBCollectionLabel<<endl;

  vector<MuonTrajectoryBuilder::TrackCand> glbTrackCands;

  event.getByToken(glbMuonsTrajToken, glbMuonsTraj);
    
  const reco::TrackCollection *glbTracks = glbMuons.product();
  
  for(unsigned int ww=0;ww<theRefits.size();ww++) {

		vector< Trajectory* > trajectories;
		vector< Trajectory* > trajectoriesUpdate;

    std::vector<std::pair<Trajectory*,reco::TrackRef> > miniMap;
    std::vector<std::pair<Trajectory*,reco::TrackRef> > miniMapUpdate;
    reco::TrackRef::key_type trackIndex = 0;
    int glbCounter = 0;


    for (reco::TrackCollection::const_iterator track = glbTracks->begin(); track!=glbTracks->end(); track++ , ++trackIndex) {

			vector< pair<Trajectory, Trajectory> > thisGlobalCombRefits; 

			reco::TrackRef glbRef(glbMuons,trackIndex);

      
			// Refit this global muon according to the refit type
			if (theRefits[ww]=="combinatoric") {
				thisGlobalCombRefits.clear();
				// Count muon hits (don't count RPC)
				int nHits = 0;
				for (trackingRecHit_iterator hit = (*track).recHitsBegin(); hit != (*track).recHitsEnd(); ++hit) {
					if (!(*hit)->isValid()) continue;
					if ((*hit)->geographicalId().det() == DetId::Tracker) continue; // skip tracker
					else if ((*hit)->geographicalId().det() == DetId::Muon) { // double check it's a muon hit
						if ((*hit)->geographicalId().subdetId() == 3) continue; // skip RPC
						nHits++;
					}
				}
				int nComb = pow(2,nHits)-1;
				for (int hp = 1; hp <= nComb; hp++) {
					pair<Trajectory, Trajectory> refits = theRefitter->refit(*track, theRefits[ww], tTopo, nHits, hp);
					if (refits.first.isValid()==false || refits.second.isValid()==false) continue;
					if (refits.first.empty()==true || refits.second.empty()==true) continue;
					thisGlobalCombRefits.push_back(refits);
				}
				//
				// Choose the best solution out of all muon hit combinatorics
				//
				if (thisGlobalCombRefits.empty()) continue;
				std::pair<Trajectory,Trajectory> bestPair;
				if (thisGlobalCombRefits.size()==1) 
					bestPair = thisGlobalCombRefits.front();
				else
					bestPair = theHighPTUtilities->select(thisGlobalCombRefits, *track);
				if (bestPair.first.empty() || bestPair.second.empty()) continue;
				
				Trajectory *refit = new Trajectory(bestPair.first);
				trajectories.push_back(refit);
				std::pair<Trajectory*,reco::TrackRef> thisPair(refit,glbRef);
				miniMap.push_back(thisPair);

				Trajectory *refitUpdate = new Trajectory(bestPair.second);
				trajectoriesUpdate.push_back(refitUpdate);
				std::pair<Trajectory*,reco::TrackRef> thisUpdatePair(refitUpdate,glbRef);
				miniMapUpdate.push_back(thisUpdatePair);
			} 
			else if (theRefits[ww]=="tracker") {
				pair<Trajectory, Trajectory> refits =theRefitter->refit(*track,theRefits[ww],tTopo,-1,-1);
				if (refits.first.isValid()==false) continue;
				if (refits.first.empty()==true) continue;

				Trajectory *refit = new Trajectory(refits.first);
				trajectories.push_back(refit);
				std::pair<Trajectory*,reco::TrackRef> thisPair(refit,glbRef);
				miniMap.push_back(thisPair);
			}
			else {
				pair<Trajectory, Trajectory> refits =theRefitter->refit(*track,theRefits[ww],tTopo,-1,-1);
				if (refits.first.isValid()==false || refits.second.isValid()==false) continue;
				if (refits.first.empty()==true || refits.second.empty()==true) continue;

				Trajectory *refit = new Trajectory(refits.first);
				trajectories.push_back(refit);
				std::pair<Trajectory*,reco::TrackRef> thisPair(refit,glbRef);
				miniMap.push_back(thisPair);

				Trajectory *refitUpdate = new Trajectory(refits.second);
				trajectoriesUpdate.push_back(refitUpdate);
				std::pair<Trajectory*,reco::TrackRef> thisUpdatePair(refitUpdate,glbRef);
				miniMapUpdate.push_back(thisUpdatePair);
			}

      if (theRefits[ww] == "dyt") dytTmp[glbCounter] = *theRefitter->getDYTInfo();
      glbCounter++;

    } // end loop on global tracks

		// Load muon-only no update tracks
		theTrackLoader->loadTracks(trajectories,event,miniMap,glbMuons, *tTopo, theRefits[ww],false);
		// Load muon-only with update tracks
		if (theRefits[ww]!="tracker") {
			theTrackLoader->loadTracks(trajectoriesUpdate,event,miniMapUpdate,glbMuons, *tTopo, theRefits[ww]+"VtxUpdate",false);
		}
		trajectories.clear();
		trajectoriesUpdate.clear();
  } // end loop on refits

  filler.insert(glbMuons, dytTmp.begin(), dytTmp.end());
  filler.fill();
  event.put(std::move(dytInfo), "dytInfo");
    
  LogTrace(metname) << "Done." << endl;    
}


