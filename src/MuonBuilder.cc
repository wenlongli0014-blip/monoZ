#include <MuonBuilder.h>

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>

#include <TRandom.h>

#include <FileInPath.h>
#include <Utils.h>


MuonBuilder::MuonBuilder(Dataset &dataset, Options const &options,
                         TabulatedRngEngine &rngEngine)
    : CollectionBuilder{dataset.Reader()},
      minPtLoose_{10.}, minPtTight_{15.},
      maxRelIsoLoose_{0.25}, maxRelIsoTight_{0.15},
      isSim_{dataset.Info().IsSimulation()},
      // Use up to 2 random numbers per muon and allow up to 5 muons before
      // repetition. This gives 10 channels for TabulatedRandomGenerator.
      tabulatedRng_{rngEngine, 10},
      srcPt_{dataset.Reader(), "Muon_pt"}, srcEta_{dataset.Reader(), "Muon_eta"},
      srcPhi_{dataset.Reader(), "Muon_phi"}, srcMass_{dataset.Reader(), "Muon_mass"},
      srcCharge_{dataset.Reader(), "Muon_charge"},
      srcIsolation_{dataset.Reader(), "Muon_pfRelIso04_all"},
      srcIsPfMuon_{dataset.Reader(), "Muon_isPFcand"},
      srcIsGlobalMuon_{dataset.Reader(), "Muon_isGlobal"},
      srcIsTrackerMuon_{dataset.Reader(), "Muon_isTracker"},
      srcIdLoose_{dataset.Reader(), "Muon_softId"},
      srcIdTight_{dataset.Reader(), "Muon_tightId"},
      srcdxy_{dataset.Reader(), "Muon_dxy"},
      srcdz_{dataset.Reader(), "Muon_dz"},
      srcTrackerLayers_{dataset.Reader(), "Muon_nTrackerLayers"} {

  if(isSim_){
    genPartId_.reset(new TTreeReaderArray<int>(dataset.Reader(), "GenPart_pdgId"));
    genPartPt_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GenPart_pt"));
    genPartEta_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GenPart_eta"));
    genPartPhi_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GenPart_phi"));
  }

  auto const roccorPath_ = Options::NodeAs<std::string>(
      options.GetConfig(), {"roccor", "path"});

  rochesterCorrection_.reset(new RoccoR(FileInPath::Resolve(roccorPath_)));
}


std::vector<Muon> const &MuonBuilder::GetLoose() const {
  Update();
  return looseMuons_;
}


std::vector<Muon> const &MuonBuilder::GetTight() const {
  Update();
  return tightMuons_;
}


void MuonBuilder::ApplyRochesterCorrection(
    int index, Muon *muon, int trackerLayers) const {

  // Apply the correction only in its domain of validity
  if (muon->p4.Pt() > 200. or std::abs(muon->p4.Eta()) > 2.4)
    return;


  double scaleFactor = 1.;

  if (isSim_) {
    auto const genMatch = FindGenMatch(*muon, 0.01);

    if (genMatch)
      scaleFactor = rochesterCorrection_->kSpreadMC(
        muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi(),  genMatch->p4.Pt());
    else
      scaleFactor = rochesterCorrection_->kSmearMC(
        muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi(),
        trackerLayers, tabulatedRng_.Rndm(2 * index));
  } else
    scaleFactor = rochesterCorrection_->kScaleDT(
      muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi());
  
  
  muon->p4.SetPtEtaPhiM(muon->p4.Pt() * scaleFactor, muon->p4.Eta(),
                        muon->p4.Phi(), muon->p4.M());
}


void MuonBuilder::Build() const {

  looseMuons_.clear();
  tightMuons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    bool const passIdLoose = srcIdLoose_[i];
    bool const passIdTight = srcIdTight_[i];

    if (not passIdLoose)
      continue;

    if (not (srcIsolation_[i] <= maxRelIsoLoose_))
      continue;

    if (not (std::abs(srcEta_[i]) < 2.4))
      continue;

    Muon muon;
    muon.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    muon.uncorrP4 = muon.p4;
    muon.charge = srcCharge_[i];

    ApplyRochesterCorrection(i, &muon, srcTrackerLayers_[i]);

    if (not (muon.p4.Pt() > minPtLoose_))
      continue;

    looseMuons_.emplace_back(muon);

    // Propagate changes in momenta of loose muons into ptmiss
    AddMomentumShift(muon.uncorrP4, muon.p4);

    if (not passIdTight)
      continue;

    if (not (srcIsolation_[i] <= maxRelIsoTight_))
      continue;

    if (not (muon.p4.Pt() > minPtTight_))
      continue;

    if (not (std::abs(srcdxy_[i]) < 0.02 and std::abs(srcdz_[i]) < 0.1))
      continue;

    tightMuons_.emplace_back(muon);
  }


  // Make sure the collections are sorted in pt
  std::sort(looseMuons_.begin(), looseMuons_.end(), PtOrdered);
  std::sort(tightMuons_.begin(), tightMuons_.end(), PtOrdered);
}


std::optional<GenParticle> MuonBuilder::FindGenMatch(
    Muon const &muon, double maxDR) const {

  unsigned iClosest = -1;
  double minDR2 = std::pow(maxDR, 2);

  for (unsigned i = 0; i < genPartId_->GetSize(); ++i) {
    if (std::abs(genPartId_->At(i)) != 13)
      // Only consider muons
      continue;

    double const dR2 = utils::DeltaR2(
      muon.p4.Eta(), muon.p4.Phi(), genPartEta_->At(i), genPartPhi_->At(i));

    if (dR2 < minDR2) {
      iClosest = i;
      minDR2 = dR2;
    }
  }

  if (iClosest != unsigned(-1)) {
    GenParticle matchedParticle{genPartId_->At(iClosest)};
    matchedParticle.p4.SetPtEtaPhiM(
      genPartPt_->At(iClosest), genPartEta_->At(iClosest), genPartPhi_->At(iClosest),
      0.1057
    );
    return matchedParticle;
  } else
    return {};
} 

