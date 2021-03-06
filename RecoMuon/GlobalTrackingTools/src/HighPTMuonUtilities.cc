#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include "RecoTracker/TransientTrackingRecHit/interface/Traj2TrackHits.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "RecoMuon/GlobalTrackingTools/interface/HighPTMuonUtilities.h"
#include "RecoMuon/TrackingTools/interface/MuonUpdatorAtVertex.h"
#include "TMath.h"

//
// Constructor
//
HighPTMuonUtilities::HighPTMuonUtilities(const edm::ParameterSet& par,
				       const MuonServiceProxy* service,
				       edm::ConsumesCollector& iC) : 
			theService(service),
			theBeamSpotInputTag(par.getParameter<edm::InputTag>("beamSpot"))
{
	edm::ParameterSet updatorPar = par.getParameter<edm::ParameterSet>("MuonUpdatorAtVertexParameters");
	theUpdatorAtVtx = new MuonUpdatorAtVertex(updatorPar,service);
  thePropagatorName = par.getParameter<std::string>("Propagator");
	theSelectorName = par.getParameter<std::string>("Selector");
	factor = par.getParameter<double>("trackRankFactor");
	curvPullCut = par.getParameter<double>("curvPullCut");
	theBeamSpotToken = iC.consumes<reco::BeamSpot>(theBeamSpotInputTag);
}

//
// Destructor
//
HighPTMuonUtilities::~HighPTMuonUtilities() {
	if (theUpdatorAtVtx) delete theUpdatorAtVtx;
}

//
// set Event
//
void HighPTMuonUtilities::setEvent(const edm::Event& event) {
  theEvent = &event;
  event.getByToken(theBeamSpotToken, beamSpot);
}

//
// Update trajectory with direction and position of global track at PCA to beam spot
//
FreeTrajectoryState 
HighPTMuonUtilities::updateMuonOnlyTraj(
		const Trajectory& refitTraj,
		const reco::Track& globalTrack) const {
	
	//std::cout << "\nHighPTMuonUtilities updateMuonOnlyTraj()\n" << std::endl;
	std::pair<bool,reco::Track> refitTrackAtPCA = 
		buildTrackFromTrajAtPCA(refitTraj, *beamSpot);
	std::pair<bool,reco::Track> globalTrackAtPCA = 
		buildTrackFromTrackAtPCA(globalTrack, *beamSpot);

	if (!refitTrackAtPCA.first || !globalTrackAtPCA.first) {
		//std::cout << "Couldn't build track from trajectory or track from track at PCA to beam spot" << std::endl;
		return FreeTrajectoryState();
	}

	std::pair<CurvilinearTrajectoryParameters,CurvilinearTrajectoryError> refitUpdateParCov = 
		KFupdateTrackWithVtx(refitTrackAtPCA.second,globalTrackAtPCA.second);
	/*
	std::cout << "par" <<std::endl;
	std::cout << refitUpdateParCov.first.vector() << std::endl;
	std::cout << "cov" <<std::endl;
	std::cout << refitUpdateParCov.second.matrix() << std::endl;
	*/

	FreeTrajectoryState updatedFTSatVtx = 
		buildFTSfromParCov(refitUpdateParCov, globalTrackAtPCA.second);

	return updatedFTSatVtx;
}

//
// Build a Track at PCA from a Trajectory
// Propagate the innermost (muon station) measurement state to surface at PCA to beam spot
// Make 6-D Helix parameters and necessary ingredients for a Track
//
//
std::pair<bool,reco::Track>
HighPTMuonUtilities::buildTrackFromTrajAtPCA(
		const Trajectory& trajectory, const reco::BeamSpot &beamSpot) const {

  const std::string metname = "Muon|RecoMuon|HighPTMuonUtilities";

  MuonPatternRecoDumper debug;
  
  TrajectoryStateOnSurface innerTSOS = trajectory.geometricalInnermostState();
	//printTSOS(innerTSOS);
  
	if (!innerTSOS.isValid()) {
		//std::cout << "\ninnerTSOS is not valid in buildTrackFromTrajAtPCA()\n" << std::endl;
		return std::make_pair(false,reco::Track());
	}

  LogTrace(metname) << "Propagate to PCA...";
	//FreeTrajectoryState ftsAtVtx = 
		//theService->propagator(thePropagatorName)->propagate(*innerTSOS.freeState(), beamSpot);  
	std::pair<bool,FreeTrajectoryState>
		extrapolationResult = theUpdatorAtVtx->propagate(innerTSOS,beamSpot);

	FreeTrajectoryState ftsAtVtx;
 if(extrapolationResult.first)
    ftsAtVtx = extrapolationResult.second;
  else{    
		//std::cout << "\nExtrapolation failed in buildTrackFromTrajAtPCA()?\n" << std::endl;
    if(TrackerBounds::isInside(innerTSOS.globalPosition())){
      //LogInfo(metname) << "Track in the Tracker: taking the innermost state instead of the state at PCA";
      ftsAtVtx = *innerTSOS.freeState();
    }
    else{
      //if ( theAllowNoVtxFlag ) {
      if (true) {
        //LogInfo(metname) << "Propagation to PCA failed, taking the innermost state instead of the state at PCA";
        ftsAtVtx = *innerTSOS.freeState();
      } else {
        //LogInfo(metname) << "Stand Alone track: this track will be rejected";
        return std::pair<bool,reco::Track>(false,reco::Track());
      }
    }
  }

    
  LogTrace(metname) << "TSOS after the extrapolation at vtx";
	//std::cout << "TSOS after the extrapolation at vtx" << std::endl;
	//std::cout << debug.dumpFTS(ftsAtVtx) << std::endl;
  
  GlobalPoint pca = extrapolationResult.second.position();
  math::XYZPoint persistentPCA(pca.x(),pca.y(),pca.z());
  GlobalVector p = extrapolationResult.second.momentum();
  math::XYZVector persistentMomentum(p.x(),p.y(),p.z());

  bool bon = true;
  if(fabs(theService->magneticField()->inTesla(GlobalPoint(0,0,0)).z()) < 0.01) bon=false;   
  double ndof = trajectory.ndof(bon);
	
	//std::cout << "HighPTMuonUtilities buildTrackFromTrajAtPCA fts q/p = " << ftsAtVtx.signedInverseMomentum() << std::endl;
	//std::cout << "HighPTMuonUtilities buildTrackFromTrajAtPCA() position = " << ftsAtVtx.position() << std::endl;
	//std::cout << "HighPTMuonUtilities buildTrackFromTrajAtPCA() momentum = " << ftsAtVtx.momentum() << std::endl;
  
	// build the track
  reco::Track track(
			trajectory.chiSquared(), 
			ndof,
			persistentPCA,
			persistentMomentum,
			ftsAtVtx.charge(),
			ftsAtVtx.curvilinearError());
	//std::cout << "HighPTMuonUtilities buildTrackFromTrajAtPCA track q/p = " << track.qoverp() << std::endl;
	
  return std::make_pair(true,track);
}

//
// Build a Transient Track at PCA from Track
//
std::pair<bool,reco::Track> 
HighPTMuonUtilities::buildTrackFromTrackAtPCA(
		const reco::Track& track, const reco::BeamSpot& beamSpot) const {
	
  // round about way to get the global track built at PCA to beam spot... 
  reco::TransientTrack transTrack_tmp(track,
				      &*theService->magneticField(),
				      theService->trackingGeometry());

  TrajectoryStateOnSurface innerGlobalTSOS = 
		transTrack_tmp.innermostMeasurementState();
	FreeTrajectoryState ftsAtVtx = 
		theService->propagator(thePropagatorName)->propagate(*innerGlobalTSOS.freeState(), beamSpot);  
  GlobalPoint pca = ftsAtVtx.position();
  math::XYZPoint persistentPCA(pca.x(),pca.y(),pca.z());

  GlobalVector p = ftsAtVtx.momentum();
  math::XYZVector persistentMomentum(p.x(),p.y(),p.z());

  reco::Track trackAtPCA(transTrack_tmp.chi2(), 
				transTrack_tmp.ndof(),
				persistentPCA,
				persistentMomentum,
				ftsAtVtx.charge(),
				ftsAtVtx.curvilinearError());

	return std::make_pair(true,trackAtPCA);
}

//
// Build a Free Trajectory State from 5D Curvilinear Parameters and Covariance 
// and x,y,z reference point from global track
//
FreeTrajectoryState
HighPTMuonUtilities::buildFTSfromParCov(
		const std::pair<CurvilinearTrajectoryParameters,CurvilinearTrajectoryError>& refitUpdateParCov,
		const reco::Track& globalTrackAtPCA) const {
	
	CurvilinearTrajectoryParameters refitUpdatePar = refitUpdateParCov.first;
	CurvilinearTrajectoryError refitUpdateCov = refitUpdateParCov.second;
	if (!refitUpdateCov.posDef()) {
		//std::cout << "Updated covariance is not positive definite" << std::endl;
		return FreeTrajectoryState();
	}
	if (refitUpdatePar.vector()[0]!=refitUpdatePar.vector()[0]) {
		//std::cout << "q/p is nan!" << std::endl;
		return FreeTrajectoryState();
	}

  GlobalPoint globalTrackPCA = GlobalPoint(
			globalTrackAtPCA.referencePoint().x(),
			globalTrackAtPCA.referencePoint().y(),
			globalTrackAtPCA.referencePoint().z());
  math::XYZPoint persistentPCA(
			globalTrackPCA.x(),globalTrackPCA.y(),globalTrackPCA.z());

	double qoverp = refitUpdatePar.vector()[0];
	double lambda = refitUpdatePar.vector()[1];
	double phi = refitUpdatePar.vector()[2];

	double pT = ROOT::Math::cos(lambda) / fabs(qoverp);
	double px = pT * ROOT::Math::cos(phi);
	double py = pT * ROOT::Math::sin(phi);
	double pz = ROOT::Math::sin(lambda) / fabs(qoverp);
  math::XYZVector persistentMomentum(px,py,pz);
	GlobalVector refitTrackP(px,py,pz);

	GlobalTrajectoryParameters refitTrajPar = GlobalTrajectoryParameters(
			globalTrackPCA, refitTrackP, refitUpdatePar.charge(), &*theService->magneticField());

	FreeTrajectoryState ftsUpdated(refitTrajPar, refitUpdateCov);
	/*
	if (ftsUpdated.hasError()) {
		std::cout << "HighPTMuonUtilities buildFTSfromParCov() ftsUpdated has Error (this is good)" << std::endl;
	}
	else {
		std::cout << "HighPTMuonUtilities buildFTSfromParCov() ftsUpdated does not have Error (this is bad)" << std::endl;
	}
	*/

	return ftsUpdated;

}

//
// Kalman Filter Update 5D Parameters and Covariance with global track
// direction and position at PCA to beam spot
//
std::pair<CurvilinearTrajectoryParameters,CurvilinearTrajectoryError>
HighPTMuonUtilities::KFupdateTrackWithVtx(
		const reco::Track &refit, const reco::Track &glbTrack) const {
	

	typedef	ROOT::Math::SMatrix<double,4,5,ROOT::Math::MatRepStd<double,4,5>> Matrix45;
	typedef	ROOT::Math::SMatrix<double,4,4,ROOT::Math::MatRepStd<double,4,4>> Matrix44;
	typedef	ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepStd<double,5,5>> Matrix55;
	typedef	ROOT::Math::SMatrix<double,5,4,ROOT::Math::MatRepStd<double,5,4>> Matrix54;
	typedef	ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5>> SymMatrix55;
	typedef	ROOT::Math::SVector<double,4> Vector4;
	typedef	ROOT::Math::SVector<double,5> Vector5;

	Matrix45 H;
	Matrix55 P;
	Matrix44 R;
	Vector4 Z_g;
	Vector5 X_r;
	Vector5 X_g;

	// Define projection matrix H
	H(0,0)=0; H(0,1)=1; H(0,2)=0; H(0,3)=0; H(0,4)=0;
	H(1,0)=0; H(1,1)=0; H(1,2)=1; H(1,3)=0; H(1,4)=0;
	H(2,0)=0; H(2,1)=0; H(2,2)=0; H(2,3)=1; H(2,4)=0;
	H(3,0)=0; H(3,1)=0; H(3,2)=0; H(3,3)=0; H(3,4)=1;
	
	/*
	std::cout << "H matrix" << std::endl;
	std::cout << H << std::endl;
	*/

	Matrix54 HT = ROOT::Math::Transpose(H);
	/*
	std::cout << "H matrix transpose" << std::endl;
	std::cout << HT << std::endl;
	*/

	Basic3DVector<double> refPt_r(refit.referencePoint());
	GlobalPoint pos_r(refPt_r);

	Basic3DVector<double> refMom_r(refit.momentum());
	GlobalVector mom_r(refMom_r);

	CurvilinearTrajectoryParameters curvil_r(
			pos_r, mom_r, refit.charge());
	//std::cout << "Refit curvilinear parameters" << std::endl;
	//std::cout << curvil_r.vector() << std::endl;

	// Get current state/cov from refit
	for (int i=0; i<5; i++) X_r(i) = curvil_r.vector()[i];
	P = refit.covariance();
	
	/*
	std::cout << "Refit prior measurement" << std::endl;
	std::cout << X_r << std::endl;
	std::cout << "Refit prior curvilinear covariance" << std::endl;
	std::cout << P << std::endl;
	*/

	Basic3DVector<double> refPt_g(glbTrack.referencePoint());
	GlobalPoint pos_g(refPt_g);

	Basic3DVector<double> refMom_g(glbTrack.momentum());
	GlobalVector mom_g(refMom_g);

	CurvilinearTrajectoryParameters curvil_g(
			pos_g, mom_g, glbTrack.charge());

	//std::cout << "Global Track curvilinear parameters" << std::endl;
	//std::cout << curvil_g.vector() << std::endl;
	for (int i=0; i<5; i++) X_g(i) = curvil_g.vector()[i];


	R = H*glbTrack.covariance()*HT;
	//std::cout << "R matrix" << std::endl;
	//std::cout << R << std::endl;

	// Calculate residuals
	Vector4 Z_r = H*X_r;
	Z_g = H*X_g;
	//std::cout << "Z_g=H*X_r" << std::endl;
	//std::cout << Z_g << std::endl;
	Vector4 y = Z_g - Z_r;
	Matrix44 S = R + H*P*HT;
	Matrix44 SInv = S;//R + H*P*HT;

	/*
	std::cout << "Z_r" << std::endl;
	std::cout << Z_r << std::endl;
	std::cout << "y=Z_g-Z_r" << std::endl;
	std::cout << y << std::endl;
	std::cout << "S" << std::endl;
	std::cout << S << std::endl;
	*/

	// Calculate Gain
	if (!SInv.Invert()) {
		std::cout << "Failed at S inversion" << std::endl;
		return std::make_pair(CurvilinearTrajectoryParameters(),CurvilinearTrajectoryError());
	}
	/*
	std::cout << "S^-1" << std::endl;
	std::cout << SInv << std::endl;
	std::cout << "S*S^-1" << std::endl;
	std::cout << S*SInv << std::endl;
	*/
	Matrix54 K = P*(HT*SInv);
	/*
	std::cout << "K" << std::endl;
	std::cout << K << std::endl;
	*/

	// Update
	AlgebraicVector5 X_u = X_r + K*y;
	Matrix55 P_u = P - K*H*P;
	
	// Print
	/*
	std::cout << "\nGlobal track parameters" << std::endl;
	std::cout << X_g << std::endl;
	std::cout << "Global track covariance" << std::endl;
	std::cout << glbTrack.covariance() << std::endl;
	std::cout << "\nMu-only track parameters before update" << std::endl;
	std::cout << X_r << std::endl;
	std::cout << "Mu-only covariance before update" << std::endl;
	std::cout << P << std::endl;
	*/
	/*
	//std::cout << "Mu-only chi2, nhits, ndof" << std::endl;
	//std::cout << refit.chi2() << " " << refit.numberOfValidHits() << " " << refit.ndof() << std::endl;
	//std::cout << "Residual" << std::endl;
	//std::cout << y << std::endl;
	//std::cout << "\nKalman gain" << std::endl;
	//std::cout << K << std::endl;
	//std::cout << "Kalman gain * residual" << std::endl;
	//std::cout << K*y << std::endl;
	//std::cout << "before cov update" << std::endl;
	//std::cout << "K" << std::endl;
	//std::cout << K << std::endl;
	//std::cout << "H" << std::endl;
	//std::cout << H << std::endl;
	//std::cout << "K*H*P" << std::endl;
	//std::cout << K*H*P << std::endl;
	//std::cout << "Covariance update" << std::endl;
	//std::cout << -K*H*P << std::endl;
	*/
	/*
	std::cout << "Updated Mu-only + vtx parameters" << std::endl;
	std::cout << X_u << std::endl;
	std::cout << "Updated Mu-only + vtx covariance" << std::endl;
	std::cout << P_u << std::endl;
	*/

	AlgebraicSymMatrix55 P_u_sym;
	for (int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			if (i>j) continue;
			P_u_sym(i,j) = P_u(i,j);
		}
	}

	CurvilinearTrajectoryParameters params(X_u);
	CurvilinearTrajectoryError cov(P_u_sym);
	/*
	// Calculate pre-fit residuals
	double prefit_vtxChi2 = ROOT::Math::Similarity(y,SInv);
	// Vertex chisquare update
	// Calculate post-fit residuals
	Vector4 pf_res = Z_g - H*X_u;
	// Calculate post-fit covariance
	Matrix44 pf_cov = R - H*P_u*HT;
	if (!pf_cov.Invert()) {
		std::cout << "Failed at post-fit covariance inversion" << std::endl;
		return std::make_pair(CurvilinearTrajectoryParameters(),CurvilinearTrajectoryError());
	}

	double vtxChi2 = ROOT::Math::Similarity(pf_res,pf_cov);
	std::cout << "\nUpdate with pos & dir at vertex" << std::endl;
	std::cout << "Mu-only + vtx update parameters" << std::endl;
	std::cout << params.vector() << std::endl;
	std::cout << "Mu-only + vtx update covariance" << std::endl;
	std::cout << cov.matrix() << std::endl;
	std::cout << "chi2 from Mu-only refit" << std::endl;
	std::cout << refit.chi2() << std::endl;
	std::cout << "pre-fit chi2 at vtx" << std::endl;
	std::cout << prefit_vtxChi2 << std::endl;
	std::cout << "\nMu-only refit + prefit vertex update chi2" << std::endl;
	std::cout << refit.chi2() + prefit_vtxChi2 << std::endl;
	std::cout << "post-fit chi2 from vtx update" << std::endl;
	std::cout << vtxChi2 << std::endl;
	std::cout << "\nMu-only refit + vertex update chi2" << std::endl;
	std::cout << refit.chi2() + vtxChi2 << std::endl;
	*/

	//std::cout << "*************" << std::endl;

	return std::make_pair(params,cov);
}

//
// Tools for selecting best trajectories among combinatoric refits
//

std::pair<Trajectory,Trajectory>
HighPTMuonUtilities::select(
		const std::vector< std::pair<Trajectory,Trajectory> >& refits,
		const reco::Track& glbTrack) const
{
	std::pair<Trajectory,Trajectory> bestPair;
	//std::cout << "Selector Name " << theSelectorName << std::endl;

	if (theSelectorName=="dxy")
		bestPair = selectBasedOnDxyPull(refits, glbTrack);
	else if (theSelectorName=="trackRank")
		bestPair = selectBasedOnTrackRank(refits,glbTrack);
	else if (theSelectorName=="curvPull")
		bestPair = selectBasedOnCurvPull(refits,glbTrack);
	else if (theSelectorName=="TEST")
		bestPair = selectBasedOnTEST(refits,glbTrack);
	else
		std::cout << "Not a valid selector" << std::endl;
	/*
	if (bestPair.first.empty() || bestPair.second.empty()) {
		std::cout << "************************************" << std::endl;
		std::cout << "Trajectories are empty after selector() this is a bad thing" << std::endl;
		std::cout << "************************************" << std::endl;
	}
	*/
	return bestPair;
}
std::pair<Trajectory,Trajectory>
HighPTMuonUtilities::selectBasedOnTrackRank(
		const std::vector< std::pair<Trajectory,Trajectory> >& refits,
		const reco::Track& glbTrack) const
{
	//std::cout << "\nChoosing trajectory based on max track rank" << std::endl;
	//std::cout << "nHits * " << factor << " - chi^2\n" << std::endl;
	int i(0);
	double bestRank = -999.;
	std::pair<Trajectory,Trajectory> bestPair;
	for (auto refit : refits) {
		Trajectory muOnlyTraj = refit.first;
		Trajectory muOnlyUpdateTraj = refit.second;
		// skip empty trajectories
		if (!muOnlyTraj.isValid() || !muOnlyUpdateTraj.isValid()) {
			i++;
			continue;
		}
		else {
			std::pair<bool,reco::Track> muOnlyTrack = 
				buildTrackFromTrajAtPCA(muOnlyTraj,*beamSpot);
			std::pair<bool,reco::Track> muOnlyUpdateTrack = 
				buildTrackFromTrajAtPCA(muOnlyUpdateTraj,*beamSpot);

			if (!muOnlyTrack.first or !muOnlyUpdateTrack.first) continue;
			double thisCurvPull = fabs((muOnlyUpdateTrack.second.qoverp() - glbTrack.qoverp())/glbTrack.qoverp());
			if (thisCurvPull > curvPullCut) continue;
			int nHits = muOnlyTraj.foundHits();
			double thisRank = nHits*factor - muOnlyUpdateTraj.chiSquared();
			//std::cout << "thisRank " << thisRank << std::endl;
			if (thisRank > bestRank) {
				bestRank = thisRank;
				bestPair = std::make_pair(muOnlyTraj,muOnlyUpdateTraj);
			}
			i++;
		}
	}
	//std::cout << "bestRank " << bestRank << std::endl;
	return bestPair;
}

//
// Smallest |(k_ref - k_glb)/k_glb|
//
std::pair<Trajectory,Trajectory>
HighPTMuonUtilities::selectBasedOnCurvPull(
		const std::vector< std::pair<Trajectory,Trajectory> >& refits,
		const reco::Track& glbTrack) const
{
	//std::cout << "\nChoosing trajectory based on smallest curvature pull" << std::endl;
	int i(0);
	double bestPull = 999.;
	std::pair<Trajectory,Trajectory> bestPair;
	for (auto refit : refits) {
		Trajectory muOnlyTraj = refit.first;
		Trajectory muOnlyUpdateTraj = refit.second;
		// skip empty trajectories
		if (!muOnlyTraj.isValid() || !muOnlyUpdateTraj.isValid()) {
			i++;
			continue;
		}
		else if (muOnlyTraj.empty() || muOnlyUpdateTraj.empty()) {
			i++;
			continue;
		}
		else {
			std::pair<bool,reco::Track> muOnlyUpdateTrack = 
				buildTrackFromTrajAtPCA(muOnlyUpdateTraj,*beamSpot);

			if (!muOnlyUpdateTrack.first) continue;
			double thisPull = fabs((muOnlyUpdateTrack.second.qoverp() - glbTrack.qoverp())/glbTrack.qoverp());
			//std::cout << "thisPull " << thisPull << std::endl;
			if (thisPull < bestPull) {
				bestPull = thisPull;
				bestPair = std::make_pair(muOnlyTraj,muOnlyUpdateTraj);
			}
			i++;
		}
	}
	//std::cout << "bestPull " << bestPull << std::endl;
	return bestPair;
}

//
// curvature pull < curvPullCut and smallest dxy
//
std::pair<Trajectory,Trajectory>
HighPTMuonUtilities::selectBasedOnDxyPull(
		const std::vector< std::pair<Trajectory,Trajectory> >& refits,
		const reco::Track& glbTrack) const
{
	//std::cout << "\nChoosing trajectory based on smallest curvature pull" << std::endl;
	int i(0);
	double bestDxyPull = 999.;
	double refDxy = glbTrack.dxy(*beamSpot);
	std::pair<Trajectory,Trajectory> bestPair;
	for (auto refit : refits) {
		Trajectory muOnlyTraj = refit.first;
		Trajectory muOnlyUpdateTraj = refit.second;
		// skip empty trajectories
		if (!muOnlyTraj.isValid() || !muOnlyUpdateTraj.isValid()) {
			i++;
			continue;
		}
		else if (muOnlyTraj.empty() || muOnlyUpdateTraj.empty()) {
			i++;
			continue;
		}
		else {
			std::pair<bool,reco::Track> muOnlyTrack = 
				buildTrackFromTrajAtPCA(muOnlyTraj,*beamSpot);
			std::pair<bool,reco::Track> muOnlyUpdateTrack = 
				buildTrackFromTrajAtPCA(muOnlyUpdateTraj,*beamSpot);

			if (!muOnlyTrack.first or !muOnlyUpdateTrack.first) continue;
			double thisCurvPull = fabs((muOnlyUpdateTrack.second.qoverp() - glbTrack.qoverp())/glbTrack.qoverp());
			if (thisCurvPull > curvPullCut) continue;
			double thisDxy = muOnlyTrack.second.dxy(*beamSpot);
			double thisDxyPull = fabs((thisDxy - refDxy)/refDxy);
			//std::cout << "thisDxyPull " << thisDxyPull << std::endl;
			if (thisDxyPull < bestDxyPull) {
				bestDxyPull = thisDxyPull;
				bestPair = std::make_pair(muOnlyTraj,muOnlyUpdateTraj);
			}
			i++;
		}
	}
	//std::cout << "bestDxyPull " << bestDxyPull << std::endl;
	return bestPair;
}

//
// Testing
//
std::pair<Trajectory,Trajectory>
HighPTMuonUtilities::selectBasedOnTEST(
		const std::vector< std::pair<Trajectory,Trajectory> >& refits,
		const reco::Track& glbTrack) const
{
	std::cout << "\nChoosing trajectory based on TEST" << std::endl;
	int i(0);
	double bestDxyPull = 999;
	double bestCurvPull = 999;
	double bestNhits = -1;
	double bestChi2 = 9999;
	double bestdptopt = 999;
	double bestprob = 0.;
	std::pair<Trajectory,Trajectory> bestPair;
	double refDxy = glbTrack.dxy(*beamSpot);
	for (auto refit : refits) {
		Trajectory muOnlyTraj = refit.first;
		Trajectory muOnlyUpdateTraj = refit.second;
		// skip empty trajectories
		if (!muOnlyTraj.isValid() || !muOnlyUpdateTraj.isValid()) {
			i++;
			continue;
		}
		else {
			std::pair<bool,reco::Track> muOnlyTrack = 
				buildTrackFromTrajAtPCA(muOnlyTraj,*beamSpot);
			std::pair<bool,reco::Track> muOnlyUpdateTrack = 
				buildTrackFromTrajAtPCA(muOnlyUpdateTraj,*beamSpot);

			double thisCurvPull = fabs((muOnlyUpdateTrack.second.qoverp() - glbTrack.qoverp())/glbTrack.qoverp());
			double thisDxy = muOnlyTrack.second.dxy(*beamSpot);
			double thisDxyPull = fabs((thisDxy - refDxy)/refDxy);
			double thisChi2 = muOnlyUpdateTrack.second.chi2();
			double thisdptopt = muOnlyUpdateTrack.second.ptError()/muOnlyUpdateTrack.second.pt();
			double thisprob = -std::log(TMath::Prob(muOnlyUpdateTrack.second.chi2(),muOnlyUpdateTrack.second.ndof()+4));
			int thisNhits = muOnlyUpdateTraj.measurements().size();

			std::cout << "a dxy pull " << thisDxyPull
				<< " a curv pull " << thisCurvPull 
				<< " a chi2 " << thisChi2
				<< " a hits " << thisNhits
				<< " a dpt/pt " << thisdptopt
				<< " a prob " << thisprob
				<< std::endl;

			if (thisCurvPull > curvPullCut) continue;
			if (thisDxyPull < bestDxyPull) {
				bestDxyPull = thisDxyPull;
				bestCurvPull = thisCurvPull;
				bestNhits = thisNhits;
				bestChi2 = thisChi2;
				bestdptopt = thisdptopt;
				bestprob = thisprob;
				bestPair = std::make_pair(muOnlyTraj,muOnlyUpdateTraj);
			}
			i++;
		}
	}
	std::cout << "chosen dxy diff " << bestDxyPull
		<< " chosen pull " << bestCurvPull
		<< " chosen chi2 " << bestChi2
		<< " chosen hits " << bestNhits
		<< " chosen dpt/pt " << bestdptopt
		<< " chosen prob " << bestprob
		<< std::endl;
	return bestPair;
}
void HighPTMuonUtilities::printTSOS(const TrajectoryStateOnSurface& tsos) const {
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
