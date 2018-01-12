#ifndef __JpsitoMuMuanalyzer_h
#define __JpsitoMuMuanalyzer_h

#include <memory>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>

// user include files                                                                                                                                           
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
const double PI=3.14159265;

using namespace std;
//histogram definition
struct histoarg{
  char name[128];		
  char title[128];
  int nbins;
  double xmin;
  double xmax;
};
enum histname{
  h_mupt,
  h_mueta,
  h_jpsimass,
  h_pv,
  h_mumdcabs,
  h_mumutrkr,
  h_mumutrkz,

  h_mumudca,
  h_mumuvtxcl,
  h_mumupt,
  h_mumumass,
  h_jpsivtxcl,
  h_jpsismass, 
Histnamesize

};

histoarg hist_arg[Histnamesize]={
  {"h_mupt","muon_pt",100,0,30},
  {"h_mueta", "Muon eta", 100, 0, 3},
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV/c^{2}]",100,1,15},
  {"h_pv","primary vertex",50,0,50},
  {"h_mumdcabs", "#mu^{-} DCA beam spot; DCA [cm]", 100, 0, 10},
  {"h_mumutrkr", "#mu^{+}#mu^{-} distance in phi-eta; [cm]", 100, 0, 50},
  {"h_mumutrkz", "#mu^{+}#mu^{-} distance in Z; [cm]", 100, 0, 100},

  {"h_mumudca",  "#mu^{+}#mu^{-} DCA; [cm]", 100, 0, 20},
  {"h_mumuvtxcl",  "#mu^{+}#mu^{-} vertex CL", 100, 0, 1},
  {"h_mumupt",    "#mu^{+}#mu^{-} pT ; pT [GeV]", 100, 0, 50},
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV/c^{2}]",
   100, 2, 20},

};

TH1F *histos[Histnamesize];
// class declaration                                                                                                                                              

class JpsitoMuMuanalyzer : public edm::EDAnalyzer{
 public:
  explicit JpsitoMuMuanalyzer(const edm::ParameterSet&);
  ~JpsitoMuMuanalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


  bool buildJpsiToMuMu(const edm::Event &);

  void calLS (double, double, double, double, double, double, double,
              double, double,  double, double, double, double, double,
              double, double, double, double, double*, double*);

  void calCosAlpha (double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double, double,
                    double, double, double, double,
		    double*, double*);
  void calCtau(RefCountedKinematicTree, double &, double &);
  double calEta(double, double, double);
  double calPhi(double, double, double);

  void clearVariables();
  bool hasBeamSpot(const edm::Event&);

  bool calClosestApproachTracks(const reco::TransientTrack,const reco::TransientTrack,double&, double&, double&);
  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
			  reco::TransientTrack &, reco::TransientTrack &,
			  double &, double &, double &, double &, double &,
			  double &, double &, double &);
  bool hasGoodJpsiMass(RefCountedKinematicTree, double &);  

  bool hasGoodJpsiVertex(const reco::TransientTrack, const reco::TransientTrack,
		       double &, double &, double &, RefCountedKinematicTree &);
  void hltReport(const edm::Event&);
  bool hasPrimaryVertex(const edm::Event &);

  void saveJpsitoMuMu(const RefCountedKinematicTree);
  void saveJpsiVertex(RefCountedKinematicTree);
  void saveJpsiCosAlpha(RefCountedKinematicTree);

  void saveJpsiLsig(RefCountedKinematicTree);
  void saveJpsiCtau(RefCountedKinematicTree);
  void saveSoftMuonVariables(pat::Muon, pat::Muon, reco::TrackRef, reco::TrackRef);
  void saveDimuVariables(double, double, double, double, double, double,double, double, double, double,double, double, double, double);
  //


  string OutputFileName_;

  bool BuildJpsiToMuMu_;

  // particle properties
  ParticleMass MuonMass_;
  float MuonMassErr_;
  double JpsiMass_;
  // labels
  
  edm::InputTag TriggerResultsLabel_;
  //  edm::EDGetTokenT<edm::TriggerResults> TriggerResultsLabel_;
  edm::InputTag BeamSpotLabel_;
  //edm::EDGetTokenT<reco::BeamSpot> BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  //edm::EDGetTokenT<reco::VertexCollection> VertexLabel_;
  edm::InputTag MuonLabel_;
  //edm::EDGetTokenT<edm::View<pat::Muon>> MuonLabel_;
  //edm::EDGetTokenT<std::vector<pat::Muon>> MuonLabel_;
  edm::InputTag TrackLabel_;
  // edm::EDGetTokenT<std::vector<pat::GenericParticle>> TrackLabel_;
  
  //edm::InputTag BeamSpotLabel_;
  // edm::InputTag VertexLabel_;
  // edm::InputTag MuonLabel_;
  vector<string> TriggerNames_;
  vector<string> LastFilterNames_;
 
 //
  double MuonMinPt_;
  double MuonMaxEta_;
  double MuonMaxDcaBs_;

  double MuMuMinPt_;
  double MuMuMinInvMass_;
  double MuMuMaxInvMass_;

  double MuMuMinVtxCl_;
  double MuMuMinLxySigmaBs_;
  double MuMuMaxDca_;
  double MuMuMinCosAlphaBs_;

  double TrkMinPt_;
  double TrkMinDcaSigBs_;
  double TrkMaxR_;
  double TrkMaxZ_;

  double JpsiMinVtxCl_;
  double JpsiMinMass_;
  double JpsiMaxMass_;
  map<string, string> mapTriggerToLastFilter_;
  reco::BeamSpot beamSpot_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  TFile* fout_;
  TTree* tree_;
  unsigned int run, event, lumiblock, nprivtx;
  vector <string> *triggernames;
  vector<int> *triggerprescales;

  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz;
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz;
  vector<double> *mumutrkr, *mumutrkz , *mumudca;
  vector<double> *mumuvtxcl, *mumulsbs, *mumulsbserr;
  vector<double> *mumucosalphabs, *mumucosalphabserr;
  vector<double> *mumumass, *mumumasserr;

  vector<bool>   *mumisgoodmuon, *mupisgoodmuon ;
  vector<int>    *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers;
  vector<int>    *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers;
  vector<double> *mumnormchi2, *mupnormchi2;
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx;
  vector<string> *mumtriglastfilter, *muptriglastfilter;
  vector<double> *mumpt, *muppt, *mumeta, *mupeta;
  vector<double> *mumtrkqual,*muptrkqual;
  int njpsi;
  vector<int>    *jpsichg; // +1 for jpsi+, -1 for bjpsi-
  vector<double> *jpsipx, *jpsipxerr, *jpsipy, *jpsipyerr, *jpsipz, *jpsipzerr, *jpsimass, *jpsimasserr;
  vector<double> *jpsivtxcl, *jpsivtxx, *jpsivtxxerr, *jpsivtxy, *jpsivtxyerr, *jpsivtxz, *jpsivtxzerr;
  vector<double> *jpsicosalphabs, *jpsicosalphabserr, *jpsilsbs, *jpsilsbserr, *jpsictau, *jpsictauerr; 
  TDatime t_begin_ , t_now_ ;
  int n_processed_, n_selected_;
  string decname;
  vector<bool> *istruemum, *istruemup, *istruetrk, *istruejpsi;

};

#endif

