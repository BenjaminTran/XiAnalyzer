//system include files
#include <memory>

#include "../interface/V0Correlation.h"

using namespace std;

V0Correlation::V0Correlation(const edm::ParameterSet& iConfig)
{
    
    //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_ass_ = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_ = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    multMax_ = iConfig.getUntrackedParameter<double>("multHigh", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multLow", 185);
    bkgnum_ = iConfig.getUntrackedParameter<int>("bkgnum", 20);
    ptbin_n_ = iConfig.getUntrackedParameter<int>("ptbin_n", 13);
    peakFactor_ = iConfig.getUntrackedParameter<double>("peakFactor", 2.0);
    sideFactor_ = iConfig.getUntrackedParameter<double>("sideFactor", 3.0);
    mis_ks_range_ = iConfig.getUntrackedParameter<double>("mis_ks_range", 0.020);
    mis_la_range_ = iConfig.getUntrackedParameter<double>("mis_la_range", 0.010);
    mis_ph_range_ = iConfig.getUntrackedParameter<double>("mis_ph_range", 0.015);

    
    sigma_ks_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_ks");
    mean_ks_ = iConfig.getUntrackedParameter<std::vector<double> >("mean_ks");
    sigma_la_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_la");
    mean_la_ = iConfig.getUntrackedParameter<std::vector<double> >("mean_la");
    ptcut_ks_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_ks");
    ptcut_la_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_la");

    rejectDaughter_ = iConfig.getUntrackedParameter<bool>("rejectDaughter");
    doRap_ = iConfig.getUntrackedParameter<bool>("doRap");
    
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _ksCollection = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("ksCollection"));
    _laCollection = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("laCollection"));
}


V0Correlation::~V0Correlation()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V0Correlation::analyze(const edm::Event& iEvent, const edm::EventSetup&
                         iSetup)
{
    using namespace edm;
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByToken(_ksCollection,v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByToken(_laCollection,v0candidates_la);
    if(!v0candidates_la.isValid()) return;
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);
    
    for(int i=0;i<ptbin_n_;i++)
    {
        pVect_trg_ks[i] = new vector<TLorentzVector>;
        pVect_trg_la[i] = new vector<TLorentzVector>;
        pVect_dau_ks[i] = new vector<TVector3>;
        pVect_dau_la[i] = new vector<TVector3>;
        pVect_trg_ks_bkg[i] = new vector<TLorentzVector>;
        pVect_trg_la_bkg[i] = new vector<TLorentzVector>;
        pVect_dau_ks_bkg[i] = new vector<TVector3>;
        pVect_dau_la_bkg[i] = new vector<TVector3>;
    }
    
    pVect_ass = new vector<TVector3>;
    

    //track selection
    int nMult_ass_good = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    hMult->Fill(nMult_ass_good);
    
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        hMult_accept->Fill(nMult_ass_good);
        for(unsigned it=0; it<v0candidates_ks->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
            
            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
            
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            double rapidity = trk.rapidity();
            double EffXchoice = 0;

            if(doRap_)
                EffXchoice = rapidity;
            else
                EffXchoice = eta;
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            //if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            //if(dlos<=5) continue;
            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);
            
            double pxd1 = dau1->px();
            double pyd1 = dau1->py();
            double pzd1 = dau1->pz();
            double pd1 = dau1->p();
            double pxd2 = dau2->px();
            double pyd2 = dau2->py();
            double pzd2 = dau2->pz();
            double pd2 = dau2->p();
            
            TVector3 dauvec1(pxd1,pyd1,pzd1);
            TVector3 dauvec2(pxd2,pyd2,pzd2);
            
            TVector3 dauvecsum(dauvec1+dauvec2);
            double v0masspiproton1 = sqrt((sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());

            double v0masspiproton2 = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))-dauvecsum.Mag2());

            //if((v0masspiproton1>=(1.115683-mis_la_range_) && v0masspiproton1<=(1.115683+mis_la_range_)) || (v0masspiproton2>=(1.115683-mis_la_range_) && v0masspiproton2<=(1.115683+mis_la_range_)) ) continue;
            
            //efficiency
            double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt));

            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();
            
            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();
            
            TLorentzVector pvector;
            pvector.SetPtEtaPhiE(pt,eta,phi,rapidity);
            
            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);
            
            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);
            
            for(int i=0;i<ptbin_n_;i++)
            {
                double rapOrEtaMaxCut = etaMax_trg_;
                double rapOrEtaMinCut = etaMin_trg_;
                double rapOrEta = eta;

                if(doRap_) 
                {
                    rapOrEtaMaxCut = 1.0;
                    rapOrEtaMinCut = -1.0;
                    rapOrEta = rapidity;
                }
                if(rapOrEta<=rapOrEtaMaxCut && rapOrEta>=rapOrEtaMinCut && pt<=ptcut_ks_[i+1] && pt>=ptcut_ks_[i]){
                    hMass_ks[i]->Fill(mass);
                    if(mass<=(mean_ks_[i]+peakFactor_*sigma_ks_[i]) && mass>=(mean_ks_[i]-peakFactor_*sigma_ks_[i])){
                        pVect_trg_ks[i]->push_back(pvector);
                        pVect_dau_ks[i]->push_back(pvector_dau1);
                        pVect_dau_ks[i]->push_back(pvector_dau2);
                        hPt_ks[i]->Fill(pt,1.0/effks);
						hEta_ks[i]->Fill(eta,1.0/effks);
                        hRap_ks[i]->Fill(rapidity);
                        hRap_ks_Lorentz[i]->Fill(pvector.E());
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_ks[i]->Fill(KET,1.0/effks);
                    }
                    if((mass<=(mean_ks_[i]-sideFactor_*sigma_ks_[i]) && mass>=0.425) || (mass<=0.57 && mass>=(mean_ks_[i]+sideFactor_*sigma_ks_[i]))){
                        pVect_trg_ks_bkg[i]->push_back(pvector);
                        pVect_dau_ks_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_ks_bkg[i]->push_back(pvector_dau2);
                        hPt_ks_bkg[i]->Fill(pt,1.0/effks);
						hEta_ks_bkg[i]->Fill(eta,1.0/effks);
                        hRap_ks_bkg[i]->Fill(rapidity);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_ks_bkg[i]->Fill(KET,1.0/effks);
                    }
                }
            }
        }
        
        for(unsigned it=0; it<v0candidates_la->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
            
            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
            
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            double rapidity = trk.rapidity();
            double EffXchoice = 0;

            if(doRap_)
                EffXchoice = rapidity;
            else
                EffXchoice = eta;
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            //if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            //if(dlos<=5) continue;
            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);

            double pxd1 = dau1->px();
            double pyd1 = dau1->py();
            double pzd1 = dau1->pz();
            double pd1 = dau1->p();
            double pxd2 = dau2->px();
            double pyd2 = dau2->py();
            double pzd2 = dau2->pz();
            double pd2 = dau2->p();
            
            TVector3 dauvec1(pxd1,pyd1,pzd1);
            TVector3 dauvec2(pxd2,pyd2,pzd2);
            
            TVector3 dauvecsum(dauvec1+dauvec2);
            double v0masspipi = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());
            double v0massee = sqrt((sqrt(0.000511*0.000511+pd1*pd1)+sqrt(0.000511*0.000511+pd2*pd2))*(sqrt(0.000511*0.000511+pd1*pd1)+sqrt(0.000511*0.000511+pd2*pd2))-dauvecsum.Mag2());
            
            //if( (v0masspipi>=(0.497614-mis_ks_range_) && v0masspipi<=(0.497614+mis_ks_range_)) || v0massee <= mis_ph_range_ ) continue;
            
            //efficiency
            double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt));

            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();
            
            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();
            
            
            TLorentzVector pvector;
            pvector.SetPtEtaPhiE(pt,eta,phi,rapidity);
            
            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);
            
            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);
            
            for(int i=0;i<ptbin_n_;i++)
            {
                double rapOrEtaMaxCut = etaMax_trg_;
                double rapOrEtaMinCut = etaMin_trg_;
                double rapOrEta = eta;

                if(doRap_)
                {
                    rapOrEtaMaxCut = 1.0;
                    rapOrEtaMinCut = -1.0;
                    rapOrEta = rapidity;
                }
                if(rapOrEta<=rapOrEtaMaxCut && rapOrEta>=rapOrEtaMinCut && pt<=ptcut_la_[i+1] && pt>=ptcut_la_[i]){
                    hMass_la[i]->Fill(mass);
                    if(mass<=(mean_la_[i]+peakFactor_*sigma_la_[i]) && mass>=(mean_la_[i]-peakFactor_*sigma_la_[i])){
                        pVect_trg_la[i]->push_back(pvector);
                        pVect_dau_la[i]->push_back(pvector_dau1);
                        pVect_dau_la[i]->push_back(pvector_dau2);
                        hPt_la[i]->Fill(pt,1.0/effla);
                        hEta_la[i]->Fill(eta,1.0/effla);
                        hRap_la[i]->Fill(rapidity);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_la[i]->Fill(KET,1.0/effla);
                    }
                    if((mass<=1.165 && mass>=(mean_la_[i]+sideFactor_*sigma_la_[i])) || (mass<=(mean_la_[i]-sideFactor_*sigma_la_[i]) && mass>=1.075)){
                        pVect_trg_la_bkg[i]->push_back(pvector);
                        pVect_dau_la_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_la_bkg[i]->push_back(pvector_dau2);
                        hPt_la_bkg[i]->Fill(pt,1.0/effla);
                        hEta_la_bkg[i]->Fill(eta,1.0/effla);
                        hRap_la_bkg[i]->Fill(rapidity);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_la_bkg[i]->Fill(KET,1.0/effla);
                    }
                }
            }
        }
        
        
        for(unsigned it=0; it<tracks->size(); ++it){
            
            const reco::Track & trk = (*tracks)[it];
            
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzvtx = trk.dz(bestvtx);
            double dxyvtx = trk.dxy(bestvtx);
            double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
            
            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError())/trk.pt()>0.10) continue;
            if(fabs(dzvtx/dzerror) > 3) continue;
            if(fabs(dxyvtx/dxyerror) > 3) continue;
            
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            if(eta<=etaMax_ass_ && eta>=etaMin_ass_ && pt<=ptMax_ass_ && pt>=ptMin_ass_) pVect_ass->push_back(pvector);
        }
        
        //Calculating signal
        int nMult_ass = (int)pVect_ass->size();
        hMult_ass->Fill(nMult_ass);
        
        for(int i=0; i<ptbin_n_; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks[i]->size();
            int nMult_trg_la = (int)pVect_trg_la[i]->size();
            
            double nMult_trg_eff_ks=0;
            double nMult_trg_eff_la=0;
            
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;
                
                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_ks = nMult_trg_eff_ks + 1.0/effks;
            }
            
			for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_la = nMult_trg_eff_la + 1.0/effla;
            }

            hMult_ks[i]->Fill(nMult_trg_ks);
            hMult_la[i]->Fill(nMult_trg_la);
            
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;
                
                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_ks[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_ks[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                    
                    if(rejectDaughter_)
                    {
                        if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }
   
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_ks/effks);///effweight_ass);
                }
            }
            
            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;
                
                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_la[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_la[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    //double pt_ass = pvector_ass.Pt();
                    
                    //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                    
		if(rejectDaughter_)
                {
                    if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                    if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                }
   
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_la/effla);///effweight_ass);
                }
            }

        }

        for(int i=0; i<ptbin_n_; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks_bkg[i]->size();
            int nMult_trg_la = (int)pVect_trg_la_bkg[i]->size();
            
            double nMult_trg_eff_ks=0;
            double nMult_trg_eff_la=0;
            
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks_bkg[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;
                
                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_ks = nMult_trg_eff_ks + 1.0/effks;
            }
            
            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la_bkg[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_la = nMult_trg_eff_la + 1.0/effla;
            }
            
            hMult_ks_bkg[i]->Fill(nMult_trg_ks);
            hMult_la_bkg[i]->Fill(nMult_trg_la);
            
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;
                
                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_ks_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_ks_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    //double pt_ass = pvector_ass.Pt();
                    
                    //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                    
                    if(rejectDaughter_)
                    {
                        if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }
   
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_ks/effks);///effweight_ass);
                }
            }
            
            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;
                
                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_la_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_la_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    //double pt_ass = pvector_ass.Pt();
                    
                    //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                   
                    if(rejectDaughter_)
                    { 
                        if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }
   
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_la/effla);///effweight_ass);
                }
            }
            
        }
        
        
        for(int i=0; i<ptbin_n_; i++)
        {
            pVectVect_trg_ks[i]->push_back(*pVect_trg_ks[i]);
            pVectVect_trg_la[i]->push_back(*pVect_trg_la[i]);
            pVectVect_dau_ks[i]->push_back(*pVect_dau_ks[i]);
            pVectVect_dau_la[i]->push_back(*pVect_dau_la[i]);
            delete pVect_trg_ks[i];
            delete pVect_trg_la[i];
            delete pVect_dau_ks[i];
            delete pVect_dau_la[i];
            pVectVect_trg_ks_bkg[i]->push_back(*pVect_trg_ks_bkg[i]);
            pVectVect_trg_la_bkg[i]->push_back(*pVect_trg_la_bkg[i]);
            pVectVect_dau_ks_bkg[i]->push_back(*pVect_dau_ks_bkg[i]);
            pVectVect_dau_la_bkg[i]->push_back(*pVect_dau_la_bkg[i]);
            delete pVect_trg_ks_bkg[i];
            delete pVect_trg_la_bkg[i];
            delete pVect_dau_ks_bkg[i];
            delete pVect_dau_la_bkg[i];
        }
        
        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);
        delete pVect_ass;
    }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
V0Correlation::beginJob()
{
    edm::Service<TFileService> fs;
    
    TH1::SetDefaultSumw2();

	/*
    edm::FileInPath fip("Demo/DemoAnalyzer/data/trkEff_pp_all_74X_origin.root");
    TFile f(fip.fullPath().c_str(),"READ");
    effhisto = (TH2F*)f.Get("rTotalEff3D");
    
	
    edm::FileInPath fip1("Demo/DemoAnalyzer/data/V0_pp13TeV_Efficiency.root");
    TFile f1(fip1.fullPath().c_str(),"READ");
    effhisto_ks = (TH2D*)f1.Get("ks_eff_1");
    effhisto_la = (TH2D*)f1.Get("la_eff_1");
	*/

    edm::FileInPath fip("XiAnalyzer/XiAnalyzer/data/Effhisto.root");
    TFile f(fip.fullPath().c_str(),"READ");
    effhisto_ks = (TH2D*)f.Get("EffHistoKs");
    effhisto_la = (TH2D*)f.Get("EffHistoLa");
    
    hMult = fs->make<TH1D>("mult",";N",600,0,600);
    hMult_ass = fs->make<TH1D>("mult_ass",";N",600,0,600);
    hMult_accept = fs->make<TH1D>("mult_acc",";N",600,0,600);
    
    for(int i=0; i<ptbin_n_; i++)
    {
        hKET_ks[i] = fs->make<TH1D>(Form("KETkshort_pt%d",i),";GeV",50000,0,25);
        hKET_la[i] = fs->make<TH1D>(Form("KETlambda_pt%d",i),";GeV",50000,0,25);
        hPt_ks[i] = fs->make<TH1D>(Form("Ptkshort_pt%d",i),";GeV",50000,0,25);
        hPt_la[i] = fs->make<TH1D>(Form("Ptlambda_pt%d",i),";GeV",50000,0,25);
        hEta_ks[i] = fs->make<TH1D>(Form("Etakshort_pt%d",i),";eta",24,-2.4,2.4);
        hRap_ks[i] = fs->make<TH1D>(Form("Rapkshort_pt%d",i),";y",30,-1.5,1.5);
        hRap_ks_Lorentz[i] = fs->make<TH1D>(Form("RapksLorentz_pt%d",i),";y",30,-1.5,1.5);
        hEta_la[i] = fs->make<TH1D>(Form("Etalambda_pt%d",i),";eta",24,-2.4,2.4);
        hRap_la[i] = fs->make<TH1D>(Form("Raplambda_pt%d",i),";y",30,-1.5,1.5);
        hMass_ks[i] = fs->make<TH1D>(Form("masskshort_pt%d",i),";GeV",2000,0,1.0);
        hMass_la[i] = fs->make<TH1D>(Form("masslambda_pt%d",i),";GeV",2000,0.5,1.5);
        hMult_ks[i] = fs->make<TH1D>(Form("mult_ks_pt%d",i),";N",250,0,250);
        hSignal_ks[i] = fs->make<TH2D>(Form("signalkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks[i] = fs->make<TH2D>(Form("backgroundkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la[i] = fs->make<TH1D>(Form("mult_la_pt%d",i),";N",250,0,250);
        hSignal_la[i] = fs->make<TH2D>(Form("signallambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la[i] = fs->make<TH2D>(Form("backgroundlambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks[i] = new vector< vector<TLorentzVector> >;
        pVectVect_trg_la[i] = new vector< vector<TLorentzVector> >;
        pVectVect_dau_ks[i] = new vector< vector<TVector3> >;
        pVectVect_dau_la[i] = new vector< vector<TVector3> >;
        hKET_ks_bkg[i] = fs->make<TH1D>(Form("KETkshort_bkg_pt%d",i),";GeV",50000,0,25);
        hKET_la_bkg[i] = fs->make<TH1D>(Form("KETlambda_bkg_pt%d",i),";GeV",50000,0,25);
        hPt_ks_bkg[i] = fs->make<TH1D>(Form("Ptkshort_bkg_pt%d",i),";GeV",50000,0,25);
        hPt_la_bkg[i] = fs->make<TH1D>(Form("Ptlambda_bkg_pt%d",i),";GeV",50000,0,25);
        hEta_ks_bkg[i] = fs->make<TH1D>(Form("Etakshort_bkg_pt%d",i),";GeV",24,-2.4,2.4);
        hRap_ks_bkg[i] = fs->make<TH1D>(Form("Rapkshort_bkg_pt%d",i),";y",30,-1.5,1.5);
        hEta_la_bkg[i] = fs->make<TH1D>(Form("Etalambda_bkg_pt%d",i),";GeV",24,-2.4,2.4);
        hRap_la_bkg[i] = fs->make<TH1D>(Form("Raplambda_bkg_pt%d",i),";y",30,-1.5,1.5);
        hMult_ks_bkg[i] = fs->make<TH1D>(Form("mult_ks_bkg_pt%d",i),";N",250,0,250);
        hSignal_ks_bkg[i] = fs->make<TH2D>(Form("signalkshort_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks_bkg[i] = fs->make<TH2D>(Form("backgroundkshort_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la_bkg[i] = fs->make<TH1D>(Form("mult_la_bkg_pt%d",i),";N",250,0,250);
        hSignal_la_bkg[i] = fs->make<TH2D>(Form("signallambda_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la_bkg[i] = fs->make<TH2D>(Form("backgroundlambda_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks_bkg[i] = new vector< vector<TLorentzVector> >;
        pVectVect_trg_la_bkg[i] = new vector< vector<TLorentzVector> >;
        pVectVect_dau_ks_bkg[i] = new vector< vector<TVector3> >;
        pVectVect_dau_la_bkg[i] = new vector< vector<TVector3> >;
    }
    
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
    
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
V0Correlation::endJob() {
    //Calculating background
    int nevttotal_ass = (int)pVectVect_ass->size();
    
    /*double etacut[7] = {-2.4,-1.6,-0.8,0,0.8,1.6,2.4};
    double kseff[6] = {1.2,1.2,1.2,1.2,1.2,1.2};
    double laeff[6] = {2,2,2,2,2,1.8};*/

    for(int i=0;i<ptbin_n_;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la[i]->size();
        
        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_ks[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_ks[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
				for(int ntrg=0;ntrg<nMult_trg;ntrg++)
				{
						TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
						double pt_trg = pvectorTmp_trg.Pt();
						double eta_trg = pvectorTmp_trg.Eta();
                        double rap_trg = pvectorTmp_trg.E();
                        double EffXchoice = 0;

                        if(doRap_)
                            EffXchoice = rap_trg;
                        else
                            EffXchoice = eta_trg;

						double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

						nMult_trg_eff = nMult_trg_eff + 1.0/effks;
				}
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;
                    
                    double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();
                    
                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        //double pt_ass = pvectorTmp_ass.Pt();
                        
                        //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                        
                        if(rejectDaughter_)
                        {
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effks);///effweight_ass);
                    }
                }
            }
        }
        
        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_la[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_la[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    nMult_trg_eff = nMult_trg_eff + 1.0/effla;
                }
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;
                    
                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();
                    
                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        //double pt_ass = pvectorTmp_ass.Pt();
                        
                        //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {                        
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effla);///effweight_ass);
                    }
                }
            }
        }

    }
    
    for(int i=0;i<ptbin_n_;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks_bkg[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la_bkg[i]->size();
        
        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_ks_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_ks_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));
                    nMult_trg_eff = nMult_trg_eff + 1.0/effks;
                }
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;
                    
                    double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();
                    
                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        //double pt_ass = pvectorTmp_ass.Pt();
                        
                        //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {                        
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effks);///effweight_ass);
                    }
                }
            }
        }
        
        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_la_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_la_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    nMult_trg_eff = nMult_trg_eff + 1.0/effla;
                }
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;
                    
                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();
                    
                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        //double pt_ass = pvectorTmp_ass.Pt();
                        
                        //double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {                        
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effla);///effweight_ass);
                    }
                }
            }
        }
        
    }
}
