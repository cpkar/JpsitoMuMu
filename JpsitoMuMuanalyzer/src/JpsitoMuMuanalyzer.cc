// -*- C++ -*-
//
// Package:    JpsiToMuMu/JpsitoMuMuanalyzer
// Class:      JpsitoMuMuanalyzer
// 
/**\class JpsitoMuMuanalyzer JpsitoMuMuanalyzer.cc JpsiToMuMu/JpsitoMuMuanalyzer/plugins/JpsitoMuMuanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Chandiprasad Kar
//         Created:  Wed, 23 Aug 2017 06:03:17 GMT
//
//


// system include files
#include "JpsitoMuMuanalyzer.h"

JpsitoMuMuanalyzer::JpsitoMuMuanalyzer(const edm::ParameterSet& iConfig):
  OutputFileName_(iConfig.getParameter<string>("OutputFileName")),
  BuildJpsiToMuMu_(iConfig.getUntrackedParameter<bool>("BuildJpsitoMuMuanalyzer")),

  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  JpsiMass_(iConfig.getUntrackedParameter<double>("JpsiMass")),

  TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
  MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),
  TriggerNames_(iConfig.getParameter< vector<string> >("TriggerNames")),
  LastFilterNames_(iConfig.getParameter< vector<string> >("LastFilterNames")),

  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),

  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")), 
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),

  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")),
  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")),

  TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
  TrkMinDcaSigBs_(iConfig.getUntrackedParameter<double>("TrkMinDcaSigBs")),
  TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
  TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),

  JpsiMinVtxCl_(iConfig.getUntrackedParameter<double>("JpsiMinVtxCl")),
  JpsiMinMass_(iConfig.getUntrackedParameter<double>("JpsiMinMass")),
  JpsiMaxMass_(iConfig.getUntrackedParameter<double>("JpsiMaxMass")),

  tree_(0),
  triggernames(0), triggerprescales(0),
  mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0),
  mupdcabs(0),  mupdcabserr(0), muppx(0),  muppy(0), muppz(0),
  mumutrkr(0), mumutrkz(0), mumudca(0),  mumuvtxcl(0),  mumulsbs(0),
  mumulsbserr(0),mumucosalphabs(0), mumucosalphabserr(0),
  mumumass(0), mumumasserr(0),
  mumisgoodmuon(0), mupisgoodmuon(0),
  mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0),
  mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0),
  mumnormchi2(0), mupnormchi2(0),mumdxyvtx(0), mupdxyvtx(0),
  mumdzvtx(0), mupdzvtx(0), mumtriglastfilter(0), muptriglastfilter(0),
  mumpt(0), muppt(0), mumeta(0), mupeta(0),mumtrkqual(0),muptrkqual(0),
  njpsi(0), jpsipx(0), jpsipxerr(0), jpsipy(0), jpsipyerr(0), jpsipz(0), jpsipzerr(0), jpsimass(0), jpsimasserr(0),
  jpsivtxcl(0), jpsivtxx(0), jpsivtxxerr(0), jpsivtxy(0), jpsivtxyerr(0), jpsivtxz(0), jpsivtxzerr(0),
  jpsicosalphabs(0), jpsicosalphabserr(0), jpsilsbs(0), jpsilsbserr(0), jpsictau(0), jpsictauerr(0)
{
   //now do what ever initialization is needed
  assert(TriggerNames_.size()==LastFilterNames_.size());
  for(size_t i=0;i<TriggerNames_.size();++i)
    mapTriggerToLastFilter_[TriggerNames_[i]]=LastFilterNames_[i];

}


JpsitoMuMuanalyzer::~JpsitoMuMuanalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
JpsitoMuMuanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   clearVariables();

   printf("\n ---------- Begin processing ---------- \n");
   n_processed_ += 1;

   run = iEvent.id().run() ;
   event = iEvent.id().event() ;
   lumiblock = iEvent.luminosityBlock();
   hltReport(iEvent);
   if ( hasBeamSpot(iEvent) ) {
     iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
     if ( bFieldHandle_.isValid() && hasPrimaryVertex(iEvent) ) {
       buildJpsiToMuMu(iEvent) ;
       n_selected_ += 1;
     }
   }
   if(njpsi>0){
     tree_->Fill();
   }
   clearVariables();
}


// ------------ method called once each job just before starting event loop  ------------
void JpsitoMuMuanalyzer::beginJob()
{
  t_begin_.Set();
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;


  fout_ = new TFile(OutputFileName_.c_str(), "RECREATE");
  fout_->cd();

  for(int i=0; i<Histnamesize; i++) {
    histos[i] = new TH1F(hist_arg[i].name, hist_arg[i].title,
			 hist_arg[i].nbins,
			 hist_arg[i].xmin, hist_arg[i].xmax);
  }
  tree_ = new TTree ("tree", "JpsiToMuMu");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i"); 
  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mumutrkr", &mumutrkr);
  tree_->Branch("mumutrkz", &mumutrkz);
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumumass", &mumumass);
  tree_->Branch("mumumasserr", &mumumasserr);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);
  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumtrkqual", &mumtrkqual); 
  tree_->Branch("muptrkqual", &muptrkqual);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("njpsi", &njpsi, "njpsi/I");
  tree_->Branch("jpsipx", &jpsipx);
  tree_->Branch("jpsipxerr", &jpsipxerr);
  tree_->Branch("jpsipy", &jpsipy);
  tree_->Branch("jpsipyerr", &jpsipyerr);
  tree_->Branch("jpsipz", &jpsipz);
  tree_->Branch("jpsipzerr", &jpsipzerr);
  tree_->Branch("jpsimass", &jpsimass);
  tree_->Branch("jpsimasserr", &jpsimasserr);
  tree_->Branch("jpsivtxcl", &jpsivtxcl);
  tree_->Branch("jpsivtxx", &jpsivtxx);
  tree_->Branch("jpsivtxxerr", &jpsivtxxerr);
  tree_->Branch("jpsivtxy", &jpsivtxy);
  tree_->Branch("jpsivtxyerr", &jpsivtxyerr);
  tree_->Branch("jpsivtxz", &jpsivtxz);
  tree_->Branch("jpsivtxzerr", &jpsivtxzerr);
  tree_->Branch("jpsicosalphabs", &jpsicosalphabs);
  tree_->Branch("jpsicosalphabserr", &jpsicosalphabserr);
  tree_->Branch("jpsilsbs", &jpsilsbs);
  tree_->Branch("jpsilsbserr", &jpsilsbserr);
  tree_->Branch("jpsictau", &jpsictau);
  tree_->Branch("jpsictauerr", &jpsictauerr);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JpsitoMuMuanalyzer::endJob() 
{
  fout_->cd();
  tree_->Write();

  for(int i = 0; i < Histnamesize; i++) {
    histos[i]->Write();
    histos[i]->Delete();
  }
  fout_->Close();

  t_now_.Set();
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();
 

}
void
JpsitoMuMuanalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// method called when ending the processing of a run  ------------
void
JpsitoMuMuanalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
void 
JpsitoMuMuanalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JpsitoMuMuanalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JpsitoMuMuanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void JpsitoMuMuanalyzer::clearVariables(){
  run=0;
  event=0;
  lumiblock=0;
  nprivtx=0;
  triggernames->clear();
  triggerprescales->clear();
  mumdcabs->clear();  mumdcabserr->clear();  mumpx->clear();   mumpy->clear();  mumpz->clear();
  mupdcabs->clear();  mupdcabserr->clear();  muppx->clear();   muppy->clear();  muppz->clear();
  mumutrkr->clear(); mumutrkz->clear();
  mumudca->clear();  mumuvtxcl->clear();   mumulsbs->clear();  mumulsbserr->clear();
  mumucosalphabs->clear();  mumucosalphabserr->clear();

  mumumass->clear();mumumasserr->clear();
  mumisgoodmuon->clear();  mupisgoodmuon->clear();
  mumnpixhits->clear();  mupnpixhits->clear();  mumnpixlayers->clear();  mupnpixlayers->clear();
  mumntrkhits->clear();  mupntrkhits->clear();  mumntrklayers->clear();  mupntrklayers->clear();

  mumnormchi2->clear(); mupnormchi2->clear();

  mumdxyvtx->clear(); mupdxyvtx->clear();
  mumdzvtx->clear(); mupdzvtx->clear();
  mumtriglastfilter->clear(); muptriglastfilter->clear();
  mumpt->clear(); muppt->clear();
  mumeta->clear(); mupeta->clear();
  mumtrkqual->clear(); muptrkqual->clear(); 
  njpsi = 0;

  jpsipx->clear(); jpsipxerr->clear(); jpsipy->clear();  jpsipyerr->clear();
  jpsipz->clear(); jpsipzerr->clear();

  jpsimass->clear(); jpsimasserr->clear();
  jpsivtxcl->clear(); jpsivtxx->clear(); jpsivtxxerr->clear(); jpsivtxy->clear(); jpsivtxyerr->clear();
  jpsivtxz->clear(); jpsivtxzerr->clear(); jpsicosalphabs->clear(); jpsicosalphabserr->clear();
  jpsilsbs->clear(); jpsilsbserr->clear(); jpsictau->clear(); jpsictauerr->clear();

}
void JpsitoMuMuanalyzer::hltReport(const edm::Event& iEvent)
{

  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try {iEvent.getByLabel( TriggerResultsLabel_, hltTriggerResults ); }
  //  try {iEvent.getByToken( TriggerResultsLabel_, hltTriggerResults ); }
  catch ( ... ) { edm::LogInfo("myHLT")
      << __LINE__ << " : couldn't get handle on HLT Trigger" ; }

  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);

    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){

      // Only consider the triggered case.                                                                                                               
      if ((*hltTriggerResults)[itrig].accept() == 1){

        string triggername = triggerNames_.triggerName(itrig);
        int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);

        // Loop over our interested HLT trigger names to find if this event contains.                                                                   
        for (unsigned int it=0; it<TriggerNames_.size(); it++){
          if (triggername.find(TriggerNames_[it]) != string::npos) {
            // save the no versioned case                                                                                                                
            triggernames->push_back(TriggerNames_[it]);
	    //cout<<triggernames<<endl;
            triggerprescales->push_back(triggerprescale);

          }}}}}

}
bool JpsitoMuMuanalyzer::hasBeamSpot(const edm::Event& iEvent){
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);

  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ;
    return false;
  }

  beamSpot_ = *beamSpotHandle;
  return true;
}
bool JpsitoMuMuanalyzer::hasPrimaryVertex(const edm::Event& iEvent){
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size();
  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) {
    primaryVertex_ = *(iVertex);
    if (primaryVertex_.isValid()) break;
  }

  if (!primaryVertex_.isValid()) return false;

  return true;
}
bool JpsitoMuMuanalyzer::buildJpsiToMuMu(const edm::Event& iEvent){
  edm::Handle< vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(MuonLabel_, patMuonHandle);
  if( patMuonHandle->size() < 2 ) return false;
  bool passed;
  reco::TransientTrack refitMupTT, refitMumTT;
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;
  double mumutrk_R, mumutrk_Z, DCAmumu;
  double MuMuLSBS, MuMuLSBSErr;
  double MuMuCosAlphaBS, MuMuCosAlphaBSErr;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  double jpsi_vtx_chisq, jpsi_vtx_cl, jpsi_mass;

  RefCountedKinematicTree vertexFitTree, barVertexFitTree;

  //--------------------
  //--------- Mu-  -----
  //--------------------
  for (vector<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin();
       iMuonM != patMuonHandle->end(); iMuonM++){

    reco::TrackRef muTrackm = iMuonM->innerTrack();
    if ( muTrackm.isNull() ) continue;

    histos[h_mupt]->Fill(muTrackm->pt());
    histos[h_mueta]->Fill(muTrackm->eta());
    if ( (muTrackm->charge() != -1) ||
	 (muTrackm->pt() < 4.0) ||
	 (fabs(muTrackm->eta()) > 2.5)) continue;

    ////// check mu- DCA to beam spot //////
    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));
    passed = hasGoodMuonDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr) ;
    histos[h_mumdcabs]->Fill(DCAmumBS);
    if ( ! passed ) continue;

    //for mu+
    for (vector<pat::Muon>::const_iterator iMuonP = patMuonHandle->begin();
	 iMuonP != patMuonHandle->end(); iMuonP++){

      reco::TrackRef muTrackp = iMuonP->innerTrack();
      if ( muTrackp.isNull() ||
	   (muTrackp->charge() != 1) ||
	   (muTrackp->pt() < 4.0) ||
	   (fabs(muTrackp->eta()) > 2.5)) continue;

      ////// check mu+ DCA to beam spot //////
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_));
      passed = hasGoodMuonDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr);
      if ( ! passed ) continue;
      if (!calClosestApproachTracks(muTrackpTT, muTrackmTT,mumutrk_R, mumutrk_Z, DCAmumu) ) continue ;
      histos[h_mumutrkr]->Fill(mumutrk_R);
      histos[h_mumutrkz]->Fill(mumutrk_Z);
      histos[h_mumudca]->Fill(DCAmumu);

      if ( mumutrk_R > TrkMaxR_ || mumutrk_Z > TrkMaxZ_ || DCAmumu > MuMuMaxDca_ ) continue;

      // check dimuon vertex /////
      passed = hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT,mu_mu_vtx_cl, mu_mu_pt,mu_mu_mass, mu_mu_mass_err,MuMuLSBS, MuMuLSBSErr,MuMuCosAlphaBS,MuMuCosAlphaBSErr);

      histos[h_mumuvtxcl]->Fill(mu_mu_vtx_cl);
      histos[h_mumupt]->Fill(mu_mu_pt);
      histos[h_mumumass]->Fill(mu_mu_mass);
      if ( !passed) continue;
      TLorentzVector mu14V, mu24V,Jpsi4V;
      mu14V.SetXYZM(iMuonM->track()->px(),iMuonM->track()->py(),iMuonM->track()->pz(),MuonMass_);
      mu24V.SetXYZM(iMuonP->track()->px(),iMuonP->track()->py(),iMuonP->track()->pz(),MuonMass_);
      Jpsi4V=mu14V+mu24V;
      if(Jpsi4V.M()<2.5 || Jpsi4V.M()>4.0) continue;
      // fit Jpsi vertex  mu- mu+ 
      if ( ! hasGoodJpsiVertex(muTrackmTT, muTrackpTT,jpsi_vtx_chisq, jpsi_vtx_cl, jpsi_mass,vertexFitTree) ) continue;
      if ( (jpsi_vtx_cl < JpsiMinVtxCl_) || (jpsi_mass > JpsiMinMass_) || (jpsi_mass < JpsiMaxMass_) ) continue;
      histos[h_jpsivtxcl]->Fill(jpsi_vtx_cl);
      histos[h_jpsismass]->Fill(jpsi_mass);
      njpsi++;
      
      saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,mumutrk_R, mumutrk_Z, DCAmumu, mu_mu_vtx_cl,MuMuLSBS, MuMuLSBSErr,MuMuCosAlphaBS, MuMuCosAlphaBSErr,mu_mu_mass, mu_mu_mass_err);
      saveSoftMuonVariables(*iMuonM, *iMuonP, muTrackm, muTrackp);
      jpsivtxcl->push_back(jpsi_vtx_cl);
      saveJpsitoMuMu(vertexFitTree);
      saveJpsiVertex(vertexFitTree);
      saveJpsiCosAlpha(vertexFitTree);
      saveJpsiLsig(vertexFitTree);
      saveJpsiCtau(vertexFitTree);

    }//mu+
  }//mu-

  if ( njpsi > 0) {
    edm::LogInfo("mu+ ") << "Found " << njpsi << "  Jpsi->mu mu ";
    return true;
  }
  return false;
}

bool  JpsitoMuMuanalyzer::hasGoodMuonDcaBs (const reco::TransientTrack muTrackTT,double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    muTrackTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot_.position().x(),beamSpot_.position().y(),beamSpot_.position().z()));
  if ( !theDCAXBS.isValid() )  return false;
  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(muDcaBs) > 2.0 )   return false;
  return true;
}

bool  JpsitoMuMuanalyzer::calClosestApproachTracks (const reco::TransientTrack trackpTT,
						    const reco::TransientTrack trackmTT,
						    double & trk_R,
						    double & trk_Z,
						    double & trk_DCA)
{
  ClosestApproachInRPhi ClosestApp;
  ClosestApp.calculate(trackpTT.initialFreeState(),
		       trackmTT.initialFreeState());
  if (! ClosestApp.status() )  return false ;
  
  GlobalPoint XingPoint = ClosestApp.crossingPoint();
  
  trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y());
  trk_Z = fabs(XingPoint.z());
  
  trk_DCA = ClosestApp.distance();
  
  return true;
}


bool JpsitoMuMuanalyzer::hasGoodMuMuVertex(const reco::TransientTrack muTrackpTT,const reco::TransientTrack muTrackmTT,
					   reco::TransientTrack &refitMupTT,
					   reco::TransientTrack &refitMumTT,
					   double & mu_mu_vtx_cl, double & mu_mu_pt,
					   double & mu_mu_mass, double & mu_mu_mass_err,
					   double & MuMuLSBS, double & MuMuLSBSErr,
					   double & MuMuCosAlphaBS,
					   double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT,
					       MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT,
					       MuonMass_,chi,ndf,MuonMassErr_));

  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);

  if ( !mumuVertexFitTree->isValid())  return false;

  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();
  if ( !mumu_KV->vertexIsValid()) return false;

  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
			     int(rint(mumu_KV->degreesOfFreedom())));

  if (mu_mu_vtx_cl < MuMuMinVtxCl_)  return false;
  mumuVertexFitTree->movePointerToTheTop();

  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();

  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();

  TLorentzVector mymum, mymup, mydimu;

  mymum.SetXYZM(refitMumTT.track().momentum().x(),
		refitMumTT.track().momentum().y(),
		refitMumTT.track().momentum().z(), MuonMass_);

  mymup.SetXYZM(refitMupTT.track().momentum().x(),
		refitMupTT.track().momentum().y(),
		refitMupTT.track().momentum().z(), MuonMass_);

  mydimu = mymum + mymup;
  mu_mu_pt = mydimu.Perp();
  mu_mu_mass = mumu_KP->currentState().mass();
  mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			matrix()(6,6));

  if ((mu_mu_pt < MuMuMinPt_) || (mu_mu_mass < MuMuMinInvMass_) || (mu_mu_mass > MuMuMaxInvMass_)) return false;
  calLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
	 beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	 mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
	 mumu_KV->error().matrix()(0,1),0.0,0.0,
	 beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	 beamSpot_.covariance()(0,1),0.0,0.0,
	 &MuMuLSBS,&MuMuLSBSErr);

  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_)  return false;

  calCosAlpha(mumu_KP->currentState().globalMomentum().x(),
	      mumu_KP->currentState().globalMomentum().y(),
	      0.0,
	      mumu_KV->position().x() - beamSpot_.position().x(),
	      mumu_KV->position().y() - beamSpot_.position().y(),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
	      mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
	      0.0,
	      0.0,
	      mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
	      mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
	      0.0,
	      mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
	      0.0,
	      0.0,
	      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);

  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_)  return false;

  return true;
}
void JpsitoMuMuanalyzer::calCosAlpha(double Vx, double Vy, double Vz,
				     double Wx, double Wy, double Wz,
				     double VxErr2, double VyErr2, double VzErr2,
				     double VxyCov, double VxzCov, double VyzCov,
				     double WxErr2, double WyErr2, double WzErr2,
				     double WxyCov, double WxzCov, double WyzCov,
				     double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( (
			  (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			  (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			  (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			  (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			  (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			  (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			 (Wnorm*Wnorm*Wnorm*Wnorm) +
			 
			 ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			  (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			  (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			  
			  (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			  (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			  (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			 (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }
}
void JpsitoMuMuanalyzer::calLS (double Vx, double Vy, double Vz,
	 double Wx, double Wy, double Wz,
	 double VxErr2, double VyErr2, double VzErr2,
	 double VxyCov, double VxzCov, double VyzCov,
	 double WxErr2, double WyErr2, double WzErr2,
	 double WxyCov, double WxzCov, double WyzCov,
	 double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


bool JpsitoMuMuanalyzer::hasGoodJpsiVertex(const reco::TransientTrack mu1TT,
					 const reco::TransientTrack mu2TT,
					 double & jpsi_vtx_chisq, double & jpsi_vtx_cl,
					 double & jpsi_mass,
					 RefCountedKinematicTree & vertexFitTree)
{

  KinematicParticleFactoryFromTransientTrack pFactory;
  float chi = 0.;
  float ndf = 0.;

  // Jpsi -> mu+ mu-
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(mu1TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(mu2TT,MuonMass_,
					      chi,ndf,MuonMassErr_));


  KinematicParticleVertexFitter fitter;
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false;
  ////cout << "particles fitted to a single vertex found: " << boolalpha << vertexFitTree->isValid() << endl; 

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex jpsi_KV = vertexFitTree->currentDecayVertex();

  if ( !jpsi_KV->vertexIsValid()) return false;
  ////cout << "Bs decay vertex found: " << boolalpha << b_KV->vertexIsValid() << endl;

  jpsi_vtx_cl = TMath::Prob((double)jpsi_KV->chiSquared(),
			 int(rint(jpsi_KV->degreesOfFreedom())));

  RefCountedKinematicParticle jpsi_KP = vertexFitTree->currentParticle();
  jpsi_mass = jpsi_KP->currentState().mass();

  //if ( (b_vtx_cl < BsMinVtxCl_) || (b_mass < BsMinMass_) || (b_mass > BsMaxMass_) ) return false;
  //printf("reco Bs cand vtxcl: %6.4f , mass: %6.4f \n", b_vtx_cl, b_mass); 

  return true;

}
void JpsitoMuMuanalyzer::saveJpsitoMuMu(const RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop(); // Bs --> phi(KK) mu+ mu-                                                                                         
  RefCountedKinematicParticle jpsi_KP = vertexFitTree->currentParticle();

  jpsipx->push_back(jpsi_KP->currentState().globalMomentum().x());
  jpsipxerr->push_back( sqrt( jpsi_KP->currentState().kinematicParametersError().matrix()(3,3) ) );
  jpsipy->push_back(jpsi_KP->currentState().globalMomentum().y());
  jpsipyerr->push_back( sqrt( jpsi_KP->currentState().kinematicParametersError().matrix()(4,4) ) );
  jpsipz->push_back(jpsi_KP->currentState().globalMomentum().z());
  jpsipzerr->push_back( sqrt( jpsi_KP->currentState().kinematicParametersError().matrix()(5,5) ) );
  jpsimass->push_back(jpsi_KP->currentState().mass());
  jpsimasserr->push_back( sqrt( jpsi_KP->currentState().kinematicParametersError().matrix()(6,6) ) );
  vertexFitTree->movePointerToTheFirstChild(); // mu1                                                                                                      
  RefCountedKinematicParticle mu1_KP = vertexFitTree->currentParticle();
  vertexFitTree->movePointerToTheNextChild();  // mu2                                                                                                         
  RefCountedKinematicParticle mu2_KP = vertexFitTree->currentParticle();
  RefCountedKinematicParticle mup_KP, mum_KP ;

  if ( mu1_KP->currentState().particleCharge() > 0 ) mup_KP = mu1_KP;
  if ( mu1_KP->currentState().particleCharge() < 0 ) mum_KP = mu1_KP;
  if ( mu2_KP->currentState().particleCharge() > 0 ) mup_KP = mu2_KP;
  if ( mu2_KP->currentState().particleCharge() < 0 ) mum_KP = mu2_KP;

  muppx->push_back(mup_KP->currentState().globalMomentum().x());
  muppy->push_back(mup_KP->currentState().globalMomentum().y());
  muppz->push_back(mup_KP->currentState().globalMomentum().z());

  mumpx->push_back(mum_KP->currentState().globalMomentum().x());
  mumpy->push_back(mum_KP->currentState().globalMomentum().y());
  mumpz->push_back(mum_KP->currentState().globalMomentum().z());

}
void
JpsitoMuMuanalyzer::saveJpsiVertex(RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex jpsi_KV = vertexFitTree->currentDecayVertex();
  jpsivtxx->push_back((*jpsi_KV).position().x());
  jpsivtxxerr->push_back(sqrt( abs(jpsi_KV->error().cxx()) ));
  jpsivtxy->push_back((*jpsi_KV).position().y());
  jpsivtxyerr->push_back(sqrt( abs(jpsi_KV->error().cyy()) ));
  jpsivtxz->push_back((*jpsi_KV).position().z());
  jpsivtxzerr->push_back(sqrt( abs(jpsi_KV->error().czz()) ));

}

void 
JpsitoMuMuanalyzer::saveJpsiCosAlpha(RefCountedKinematicTree vertexFitTree)
{
  // alpha is the angle in the transverse plane between the B0 momentum                                                                             
  // and the seperation between the B0 vertex and the beamspot                                                                                           

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle jpsi_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   jpsi_KV = vertexFitTree->currentDecayVertex();

  double cosAlphaBS, cosAlphaBSErr;

  calCosAlpha(jpsi_KP->currentState().globalMomentum().x(),
              jpsi_KP->currentState().globalMomentum().y(),
              jpsi_KP->currentState().globalMomentum().z(),
              jpsi_KV->position().x() - beamSpot_.position().x(),
              jpsi_KV->position().y() - beamSpot_.position().y(),
              jpsi_KV->position().z() - beamSpot_.position().z(),
              jpsi_KP->currentState().kinematicParametersError().matrix()(3,3),
              jpsi_KP->currentState().kinematicParametersError().matrix()(4,4),
              jpsi_KP->currentState().kinematicParametersError().matrix()(5,5),
              jpsi_KP->currentState().kinematicParametersError().matrix()(3,4),
              jpsi_KP->currentState().kinematicParametersError().matrix()(3,5),
              jpsi_KP->currentState().kinematicParametersError().matrix()(4,5),
              jpsi_KV->error().cxx() + beamSpot_.covariance()(0,0),
              jpsi_KV->error().cyy() + beamSpot_.covariance()(1,1),
              jpsi_KV->error().czz() + beamSpot_.covariance()(2,2),
              jpsi_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
              jpsi_KV->error().matrix()(0,2) + beamSpot_.covariance()(0,2),
              jpsi_KV->error().matrix()(1,2) + beamSpot_.covariance()(1,2),
              &cosAlphaBS,&cosAlphaBSErr);

  jpsicosalphabs->push_back(cosAlphaBS);
  jpsicosalphabserr->push_back(cosAlphaBSErr);

}
void 
JpsitoMuMuanalyzer::saveJpsiLsig(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex jpsi_KV = vertexFitTree->currentDecayVertex();
  double LSBS, LSBSErr;

  calLS (jpsi_KV->position().x(), jpsi_KV->position().y(), 0.0,
         beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
         jpsi_KV->error().cxx(), jpsi_KV->error().cyy(), 0.0,
         jpsi_KV->error().matrix()(0,1), 0.0, 0.0,
         beamSpot_.covariance()(0,0), beamSpot_.covariance()(1,1), 0.0,
         beamSpot_.covariance()(0,1), 0.0, 0.0,
         &LSBS,&LSBSErr);

  jpsilsbs->push_back(LSBS);
  jpsilsbserr->push_back(LSBSErr);

}

void
JpsitoMuMuanalyzer::calCtau(RefCountedKinematicTree vertexFitTree,
		     double &jpsictau, double &jpsictauerr)
{
  //calculate ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)                                                                                                     

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle jpsi_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   jpsi_KV = vertexFitTree->currentDecayVertex();

  double betagamma = (jpsi_KP->currentState().globalMomentum().mag()/JpsiMass_);

  // calculate ctau error. Momentum error is negligible compared to                                                                                       
  // the vertex errors, so don't worry about it                                                                                                         

  GlobalPoint BVP = GlobalPoint( jpsi_KV->position() );
  GlobalPoint PVP = GlobalPoint( primaryVertex_.position().x(),
                                 primaryVertex_.position().y(),
                                 primaryVertex_.position().z() );
  GlobalVector sep3D = BVP-PVP;
  GlobalVector pBV = jpsi_KP->currentState().globalMomentum();
  jpsictau = (JpsiMass_* (sep3D.dot(pBV)))/(pBV.dot(pBV));

  GlobalError BVE = jpsi_KV->error();
  GlobalError PVE = GlobalError( primaryVertex_.error() );
  VertexDistance3D theVertexDistance3D;
  Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVP, PVE) );
  double myError = TheMeasurement.error();
  double scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );
  jpsictauerr =  (myError*scale)/betagamma;

}
void 
JpsitoMuMuanalyzer::saveJpsiCtau(RefCountedKinematicTree vertexFitTree)
{
  double jpsictau_temp, jpsictauerr_temp;
  calCtau(vertexFitTree, jpsictau_temp, jpsictauerr_temp);
  jpsictau->push_back(jpsictau_temp);
  jpsictauerr->push_back(jpsictauerr_temp);
}
double
JpsitoMuMuanalyzer::calEta (double Px, double Py, double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double
JpsitoMuMuanalyzer::calPhi (double Px, double Py, double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

void
JpsitoMuMuanalyzer::saveSoftMuonVariables(pat::Muon iMuonM, pat::Muon iMuonP,
				   reco::TrackRef muTrackm, reco::TrackRef muTrackp)
{

  mumisgoodmuon->push_back(muon::isGoodMuon(iMuonM, muon::TMOneStationTight));
  mupisgoodmuon->push_back(muon::isGoodMuon(iMuonP, muon::TMOneStationTight));
  mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
  mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
  mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());
  mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());

  mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
  mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
  mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
  mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());

  mumnormchi2->push_back(muTrackm->normalizedChi2());
  mupnormchi2->push_back(muTrackp->normalizedChi2());

  mumtrkqual->push_back(muTrackm->quality(reco::TrackBase::highPurity));   /* added */
  muptrkqual->push_back(muTrackp->quality(reco::TrackBase::highPurity));

  mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position()));
  mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position()));

  mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position()));
  mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position()));

  mumpt->push_back(muTrackm->pt());
  muppt->push_back(muTrackp->pt());

  mumeta->push_back(muTrackm->eta());
  mupeta->push_back(muTrackp->eta());

}
void
JpsitoMuMuanalyzer::saveDimuVariables(double DCAmumBS, double DCAmumBSErr,
			       double DCAmupBS, double DCAmupBSErr,
			       double mumutrk_R, double mumutrk_Z,
			       double DCAmumu,  double mu_mu_vtx_cl,
			       double MuMuLSBS, double MuMuLSBSErr,
			       double MuMuCosAlphaBS, double MuMuCosAlphaBSErr,
			       double mu_mu_mass, double mu_mu_mass_err)

{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);

  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);

  mumutrkr->push_back(mumutrk_R);
  mumutrkz->push_back(mumutrk_Z);
  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl);
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr);

  mumumass->push_back(mu_mu_mass);
  mumumasserr->push_back(mu_mu_mass_err);
}



//define this as a plug-in
DEFINE_FWK_MODULE(JpsitoMuMuanalyzer);
