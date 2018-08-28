#ifndef RecoMuon_GlobalMuonProducer_HighPTMuonProducer_H
#define RecoMuon_GlobalMuonProducer_HighPTMuonProducer_H

/**  \class HighPTMuonProducer
 * 
 *   High PT muon reconstructor:
 *   Refits global tracks using DT, CSC, RPC information
 *
 *
 *
 *
 *
 *   \author	Christian Schnaible (UCLA)
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonRefitter.h"
#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonUtilities.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackLoader.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
// Input and output collection

#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

typedef edm::ValueMap<reco::DYTInfo> DYTestimators;

namespace edm {class ParameterSet; class Event; class EventSetup;}

class MuonTrackFinder;
class MuonServiceProxy;

class HighPTMuonProducer : public edm::stream::EDProducer<> {

	public:

		/// constructor with config
		HighPTMuonProducer(const edm::ParameterSet&);
		
		/// destructor
		virtual ~HighPTMuonProducer(); 
		
		/// reconstruct muons
		virtual void produce(edm::Event&, const edm::EventSetup&) override;
		
	private:

	/// STA Label
	edm::InputTag theGLBCollectionLabel;
	edm::InputTag theRECOMuonCollectionLabel;
	edm::InputTag theGENParticleCollectionLabel;
	edm::EDGetTokenT<reco::TrackCollection> glbMuonsToken;
	edm::EDGetTokenT< std::vector<reco::Muon> > recoMuonsToken;
	edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken;
	edm::InputTag theBeamSpotInputTag;
	edm::EDGetTokenT< reco::BeamSpot > theBeamSpotToken;
	double dRcut;
	edm::EDGetTokenT<std::vector<Trajectory> > glbMuonsTrajToken;

	/// the event setup proxy, it takes care the services update
	MuonServiceProxy* theService;
	HighPTMuonUtilities* theHighPTUtilities;

	HighPTMuonRefitter* theRefitter;

	MuonTrackLoader* theTrackLoader;

	std::string theAlias;
	std::vector<std::string> theRefits;
	std::string theSelectorName;
	std::string BaseTrackType;

	void setAlias( std::string alias ){
		alias.erase( alias.size() - 1, alias.size() );
		theAlias=alias;
	}
  
};

#endif
