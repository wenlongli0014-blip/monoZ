#include <ZGammaTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TVector2.h>

#include <HZZException.h>
#include <Utils.h>



namespace po = boost::program_options;


int const ZGammaTrees::maxSize_;


ZGammaTrees::ZGammaTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      srcRun_{dataset.Reader(), "run"},
      srcLumi_{dataset.Reader(), "luminosityBlock"},
      srcEvent_{dataset.Reader(), "event"},
      photonBuilder_{dataset},
      photonPrescales_{dataset, options},
      photonWeight_{dataset, options, &photonBuilder_},
      triggerFilter_{dataset, options, &runSampler_},
      gJetsWeight_{dataset, &photonBuilder_},
      //photonFilter_{dataset, options},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  if (isSim_) {
    srcLHEVpt_.reset(new TTreeReaderValue<Float_t>(dataset.Reader(), "LHE_Vpt"));

    numGenPart_.reset(new TTreeReaderValue<UInt_t>(dataset.Reader(), "nGenPart"));
    genPartPdgId_.reset(new TTreeReaderArray<Int_t>(dataset.Reader(), "GenPart_pdgId"));
    genPartPt_.reset(new TTreeReaderArray<Float_t>(dataset.Reader(), "GenPart_pt"));
    genPartEta_.reset(new TTreeReaderArray<Float_t>(dataset.Reader(), "GenPart_eta"));
    genPartPhi_.reset(new TTreeReaderArray<Float_t>(dataset.Reader(), "GenPart_phi"));
    genPartStatus_.reset(new TTreeReaderArray<Int_t>(dataset.Reader(), "GenPart_status"));
    genPartStatusFlags_.reset(new TTreeReaderArray<Int_t>(dataset.Reader(), "GenPart_statusFlags"));
  }

  tauBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_, &photonBuilder_});
  ptMissBuilder_.PullCalibration({&photonBuilder_});

  weightCollector_.Add(&photonWeight_);
  weightCollector_.Add(&gJetsWeight_);

  CreateWeightBranches();

  std::cout<<dataset.Info().Name()<<std::endl;
  datasetName_ = dataset.Info().Name();
  AddBranch("lepton_cat", &leptonCat_);
  AddBranch("jet_cat", &jetCat_);
  AddBranch("jet_size", &jetSize_);
  AddBranch("photon_pt", &photonPt_);
  AddBranch("photon_eta", &photonEta_);
  AddBranch("photon_phi", &photonPhi_);
  AddBranch("photon_mass", &photonMass_);
  AddBranch("photon_r9", &photonR9_);
  AddBranch("photon_sieie", &photonSieie_);
  AddBranch("ptmiss", &missPt_);
  AddBranch("ptmiss_phi", &missPhi_);
  AddBranch("mT", &mT_);

  AddBranch("lepton_pt", leptonPt_, "lepton_pt[2]/F");
  AddBranch("lepton_eta", leptonEta_, "lepton_eta[2]/F");
  AddBranch("lepton_phi", leptonPhi_, "lepton_phi[2]/F");

  AddBranch("ll_pt", &llPt_);
  AddBranch("num_pv_good", &numPVGood_);
  AddBranch("trigger_weight", &triggerWeight_);
  // AddBranch("photon_reweighting", &photonReweighting_);
  // AddBranch("photon_nvtx_reweighting", &photonNvtxReweighting_);
  // AddBranch("photon_eta_reweighting", &photonEtaReweighting_);
  // AddBranch("mean_weight", &meanWeight_);
  // AddBranch("is_overlapped", &isOverlapped_);
  // AddBranch("sm_DjjVBF", &smDjjVBF_);
  // AddBranch("a2_DjjVBF", &a2DjjVBF_);
  // AddBranch("a3_DjjVBF", &a3DjjVBF_);
  // AddBranch("L1_DjjVBF", &l1DjjVBF_);

  if (storeMoreVariables_) {
    AddBranch("run", &run_);
    AddBranch("lumi", &lumi_);
    AddBranch("event", &event_);
  }

  datasetLHEVptUpperLimitInc_.reset();
  auto const LHEVptUpperLimitIncSettingsNode = dataset.Info().Parameters()["LHE_Vpt_upper_limit_inc"];
  if (LHEVptUpperLimitIncSettingsNode and not LHEVptUpperLimitIncSettingsNode.IsNull()) {
    datasetLHEVptUpperLimitInc_.emplace(LHEVptUpperLimitIncSettingsNode.as<Float_t>());
  }

  auto const WGSettingsNode = dataset.Info().Parameters()["wgamma_lnugamma"];
  if (WGSettingsNode and not WGSettingsNode.IsNull()) {
    labelWGamma_ = WGSettingsNode.as<std::string>();
    genPhotonBuilder_.emplace(dataset);
  }
  auto const ZGSettingsNode = dataset.Info().Parameters()["zgamma_2nugamma"];
  if (ZGSettingsNode and not ZGSettingsNode.IsNull()) {
    labelZGamma_ = ZGSettingsNode.as<std::string>();
    genPhotonBuilder_.emplace(dataset);
  }

  auto const &isQCDNode = dataset.Info().Parameters()["mc_qcd"];
  isQCD_ = (isQCDNode and not isQCDNode.IsNull() and isQCDNode.as<bool>());

  isZGToLLG_ = (datasetName_.rfind("ZGTo2LG", 0) == 0 || datasetName_.rfind("ZGToLLG", 0) == 0);
  isDYJetsToLL_ = (datasetName_.rfind("DYJetsToLL", 0) == 0);

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  v_jetCat_ = {"_eq0jets","_eq1jets","_geq2jets"};
  v_analysisCat_ = {"_eq0jets","_eq1jets","_geq2jets_discrbin1",
    "_geq2jets_discrbin2","_geq2jets_discrbin3","_geq2jets_discrbin4",
    "_geq2jets_discrbin5","_geq2jets_discrbin6","_geq2jets_discrbin7"};
  applyNvtxWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_nvtx_reweighting"});
  applyEtaWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_eta_reweighting"});
  applyPtWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_pt_reweighting"});
  applyMassLineshape_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_mass_lineshape"});
  applyMeanWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_mean_weights"});
  utils::loadInstrMETWeights(applyNvtxWeights_, applyEtaWeights_, applyPtWeights_, applyMassLineshape_, nVtxWeight_map_, etaWeight_map_, ptWeight_map_, lineshapeMassWeight_map_, v_jetCat_, options);
  utils::loadMeanWeights(applyMeanWeights_, meanWeight_map_, v_analysisCat_, options);
}


po::options_description ZGammaTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool ZGammaTrees::ProcessEvent() {
  if (not ApplyCommonFilters())
    return false;

  auto const photon = CheckPhotons();
  if (photon == nullptr)
    return false;
  //if (datasetName_ == "SinglePhoton" && !isSim_)
    //if (triggerFilter_("photon"))
      //return false;
  auto const leptonResult = CheckLeptons();
  if (not leptonResult)
    return false;

  auto const &[leptonCat, l1, l2] = leptonResult.value();
  if(datasetName_ != "SinglePhoton") {
    switch (leptonCat) {
      case LeptonCat::kEE:
        if (not triggerFilter_("ee"))
          return false;
        break;
      case LeptonCat::kMuMu:
        if (not triggerFilter_("mumu"))
          return false;
        break;
    }
  }
  else if(datasetName_ == "SinglePhoton") {
    if (triggerFilter_("ee") || triggerFilter_("mumu"))
      return false;
  }
  leptonCat_ = int(leptonCat);
  TLorentzVector const p4LL = l1->p4 + l2->p4;
  llPt_ = p4LL.Pt();
  llEta_ = p4LL.Eta();
  llPhi_ = p4LL.Phi();
  llMass_ = p4LL.M();

  if (std::abs(p4LL.M() - kNominalMZ_) > zMassWindow_)
    return false;

  if (datasetLHEVptUpperLimitInc_.has_value() and not (*srcLHEVpt_->Get() <= datasetLHEVptUpperLimitInc_.value()))
    return false;

  // Avoid double counting for ZGamma overlap between 2 samples
  if (labelZGamma_ != "" and genPhotonBuilder_) {
    double genPhotonPt = genPhotonBuilder_->P4Gamma().Pt();
    if (labelZGamma_ == "inclusive" and genPhotonPt >= 130)
      return false;
    if (labelZGamma_ == "pt130" and genPhotonPt < 130)
      return false;
  }

  // Avoid double counting for WGamma overlap between 2 samples
  if (labelWGamma_ != "" and genPhotonBuilder_) {
    double genPhotonPt = genPhotonBuilder_->P4Gamma().Pt();
    if (labelWGamma_ == "inclusive" and genPhotonPt >= 130)
      return false;
    if (labelWGamma_ == "pt130" and (genPhotonPt < 130 or genPhotonPt >= 500))
      return false;
    if (labelWGamma_ == "pt500" and genPhotonPt < 500)
      return false;
  }

  // Resolve G+jet/QCD mixing (avoid double counting of photons):
  // QCD samples allow prompt photons of pT > 10, for gamma+jets it's 25
  if (isQCD_) {
    auto const &photons = photonBuilder_.Get();
    for (int i = 0 ; i < int(photons.size()) ; i++) {
      if (photons[i].flavour == Photon::Origin::PromptPhoton 
            and photons[i].genP4.Pt() > 25.) {
        return false;  // Remove all QCD events with a photon of pT > 25
      }
    }
  }

  if (photon->p4.Pt() < minPtLL_)
    return false;


  isOverlapped_ = false;

  if (isSim_ && (isZGToLLG_ || isDYJetsToLL_)) {
    for (unsigned i = 0; i < genPartPdgId_->GetSize(); ++i) {
      // Particle is in the final state
      if (not (genPartStatus_->At(i) == 1))
        continue;

      // is a photon
      if (not (genPartPdgId_->At(i) == 22))
        continue;

      // isPrompt or fromHardProcess
      if (not ((genPartStatusFlags_->At(i) & 1) || (genPartStatusFlags_->At(i) >> 8) & 1))
        continue;

      if (not (genPartPt_->At(i) > 15 && std::abs(genPartEta_->At(i) < 2.6)))
        continue;

      // std::cout << "Gen photon!" << std::endl;

      bool gen_photon_isolated = true;

      for (unsigned j = 0; j < genPartPdgId_->GetSize(); ++j) {
        // Other particles
        if (j == i)
          continue;

        // is in the final state
        if (not (genPartStatus_->At(j) == 1))
          continue;

        // is NOT a photon
        if (genPartPdgId_->At(j) == 22)
          continue;

        // fromHardProcess
        if (not ((genPartStatusFlags_->At(j) >> 8) & 1))
          continue;

        if (not (genPartPt_->At(j) > 5))
          continue;

        // std::cout << std::sqrt(utils::DeltaR2(genPartEta_->At(i), genPartPhi_->At(i), genPartEta_->At(j), genPartPhi_->At(j))) << std::endl;

        if (std::sqrt(utils::DeltaR2(genPartEta_->At(i), genPartPhi_->At(i), genPartEta_->At(j), genPartPhi_->At(j))) < 0.05) {
          gen_photon_isolated = false;
          break;
        }
      }

      if (gen_photon_isolated) {
        // std::cout << "Event removed!" << std::endl;
        isOverlapped_ = true;
      }
    }
  }

  if (isSim_ && isDYJetsToLL_) {
    if (isOverlapped_) {
      return false;
    }
  }


  auto const p4Miss = ptMissBuilder_.Get().p4 + p4LL;
  missPt_ = p4Miss.Pt();
  missPhi_ = p4Miss.Phi();

  if (std::abs(TVector2::Phi_mpi_pi(photon->p4.Phi() - p4Miss.Phi())) 
        < minDphiLLPtMiss_)
    return false;


  auto const &jets = jetBuilder_.Get();

  for (auto const &jet : jets) {
    if (bTagger_(jet))
      return false;

    if (std::abs(TVector2::Phi_mpi_pi(jet.p4.Phi() - p4Miss.Phi())) 
          < minDphiJetPtMiss_)
      return false;
  }

  // if (DPhiPtMiss2(p4Miss, {&jetBuilder_, &photonBuilder_}) < minDphiLeptonsJetsPtMiss_)
  //   return false;

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (jets.size() == 1)
    jetCat_ = int(JetCat::kEq1J);
  else
    jetCat_ = int(JetCat::kGEq2J);

  analysisCat_ = jetCat_;

  // Veto event if the photon fails the following identification
  if (!photon->passElecVeto)
    return false;

  if (photon->sieie < 0.001)
    return false;

  photonSieie_ = photon->sieie;
  photonR9_ = photon->r9;
  // // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  // // Reweighting
  // photonReweighting_ = 1.;
  // photonNvtxReweighting_ = 1.;
  // photonEtaReweighting_ = 1.;
  // // In nvtx
  // if (applyNvtxWeights_) {
  //   std::map<double, std::pair<double,double> >::iterator itlow;
  //   itlow = nVtxWeight_map_["_ll"].upper_bound(*srcNumPVGood_);
  //   if (itlow == nVtxWeight_map_["_ll"].begin()) 
  //     throw HZZException{
  //       "You are trying to access your NVtx reweighting map outside of bin "
  //       "boundaries."
  //     };
  //   itlow--;
  //   photonReweighting_ *= itlow->second.first;
  //   photonNvtxReweighting_ *= itlow->second.first;
  // }
  // // In eta
  // if (applyEtaWeights_ and jetCat_ == int(JetCat::kGEq2J)) { // Don't apply it for 0 and 1 jet categories
  //   std::map<double, std::pair<double,double> >::iterator itlow;
  //   itlow = etaWeight_map_["_ll"+v_jetCat_[jetCat_]].upper_bound(fabs(photon->p4.Eta())); //look at which bin in the map currentEvt.eta corresponds
  //   if (itlow == etaWeight_map_["_ll" + v_jetCat_[jetCat_]].begin())
  //     throw HZZException{
  //       "You are trying to access your Eta reweighting map outside of bin "
  //       "boundaries."
  //     };
  //   itlow--;
  //   photonReweighting_ *= itlow->second.first;
  //   photonEtaReweighting_ *= itlow->second.first;
  // }
  // // In pT
  // if (applyPtWeights_) {
  //   std::map<double, std::pair<double,double> >::iterator itlow;
  //   itlow = ptWeight_map_["_ll"+v_jetCat_[jetCat_]].upper_bound(photon->p4.Pt()); //look at which bin in the map currentEvt.pT corresponds
  //   if (itlow == ptWeight_map_["_ll" + v_jetCat_[jetCat_]].begin())
  //     throw HZZException{
  //       "You are trying to access your Pt reweighting map outside of bin "
  //       "boundaries."
  //     };
  //   itlow--;
  //   photonReweighting_ *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
  // }

  // Give mass to photon
  double photonMass = 0;
  if (applyMassLineshape_) {
    photonMass = lineshapeMassWeight_map_["_ll"]->GetRandom();
  }
  TLorentzVector photonWithMass;
  photonWithMass.SetPtEtaPhiM(photon->p4.Pt(), photon->p4.Eta(), 
    photon->p4.Phi(), photonMass);

  photonPt_ = photonWithMass.Pt();
  photonEta_ = photonWithMass.Eta();
  photonPhi_ = photonWithMass.Phi();
  photonMass_ = photonWithMass.M();


  double const eT =
      std::sqrt(std::pow(photonWithMass.Pt(), 2)
          + std::pow(photonWithMass.M(), 2))
      + std::sqrt(std::pow(p4Miss.Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(
      std::pow(eT, 2) - std::pow((photonWithMass + p4Miss).Pt(), 2));

  triggerWeight_ = 1.0;
  if (not isSim_) {
    if(datasetName_ == "SinglePhoton")
      triggerWeight_ = photonPrescales_.GetPhotonPrescale(photonWithMass.Pt());
    if (triggerWeight_ == 0)
      return false;
  }

  // Apply the event veto for photons failing an OR of the addtional IDs of photon PF ID, sigmaIPhiIPhi > 0.001,
  //  MIPTotalEnergy < 4.9, |seedtime| < 2ns, and seedtime < 1ns for 2018
  //if (!photonFilter_())
    //return false;

  numPVGood_ = *srcNumPVGood_;

  // auto const &djjVBF = vbfDiscriminant_.Get(photonWithMass, p4Miss, jets);
  // smDjjVBF_ = djjVBF.at(VBFDiscriminant::DjjVBF::SM);
  // a2DjjVBF_ = djjVBF.at(VBFDiscriminant::DjjVBF::a2);
  // a3DjjVBF_ = djjVBF.at(VBFDiscriminant::DjjVBF::a3);
  // l1DjjVBF_ = djjVBF.at(VBFDiscriminant::DjjVBF::L1);

  // if (jetCat_ == int(JetCat::kGEq2J)) {
  //   if (smDjjVBF_ >= 0. and smDjjVBF_ < 0.05) analysisCat_ = 2;
  //   else if (smDjjVBF_ >= 0.05 and smDjjVBF_ < 0.1) analysisCat_ = 3;
  //   else if (smDjjVBF_ >= 0.1 and smDjjVBF_ < 0.2) analysisCat_ = 4;
  //   else if (smDjjVBF_ >= 0.2 and smDjjVBF_ < 0.8) analysisCat_ = 5;
  //   else if (smDjjVBF_ >= 0.8 and smDjjVBF_ < 0.9) analysisCat_ = 6;
  //   else if (smDjjVBF_ >= 0.9 and smDjjVBF_ < 0.95) analysisCat_ = 7;
  //   else if (smDjjVBF_ >= 0.95) analysisCat_ = 8;
  // }

  // // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  // // Get mean weights
  // if (applyMeanWeights_) {
  //   std::map<double, double>::iterator itlow;
  //   itlow = meanWeight_map_[v_analysisCat_[analysisCat_]].upper_bound(mT_); //look at which bin in the map mt corresponds
  //   if (itlow == meanWeight_map_[v_analysisCat_[analysisCat_]].begin())
  //     throw HZZException{
  //       "You are trying to access your mean weight map outside of bin "
  //       "boundaries."
  //     };
  //   itlow--;
  //   meanWeight_ = itlow->second;
  // }

  // if (applyMeanWeights_ and meanWeight_ == 0)
  //   return false;

  //if (storeMoreVariables_)
  //FillMoreVariables(jets);

  std::array<const Lepton *, 2> leptons = {l1, l2};

  for (int i = 0; i < 2; ++i) {
    leptonPt_[i] = leptons[i]->p4.Pt();
    leptonEta_[i] = leptons[i]->p4.Eta();
    leptonPhi_[i] = leptons[i]->p4.Phi();
  }

  jetSize_ = jets.size();
  FillTree();
  return true;

}


Photon const *ZGammaTrees::CheckPhotons() const {
  auto const &photons = photonBuilder_.Get();

  if (photons.size() != 1)
    return nullptr;
  if (tauBuilder_.Get().size() > 0)
    return nullptr;

  return &photons[0];
}

std::optional<std::tuple<ZGammaTrees::LeptonCat, Lepton const *, Lepton const *>>
ZGammaTrees::CheckLeptons() const {
  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  bool tag = (looseElectrons.size() == 0 && looseMuons.size() == 2) || 
      (looseElectrons.size() == 2 && looseMuons.size() == 0);
  if (!tag)
    return {};

  LeptonCat leptonCat;
  Lepton const *l1, *l2;

  if (tightElectrons.size() == 2) {
    leptonCat = LeptonCat::kEE;
    l1 = &tightElectrons[0];
    l2 = &tightElectrons[1];
  } else if (tightMuons.size() == 2) {
    leptonCat = LeptonCat::kMuMu;
    l1 = &tightMuons[0];
    l2 = &tightMuons[1];
  } else
    return {};

  if (l1->charge * l2->charge > 0)
    return {};

  if (not (l1->p4.Pt() > 15. && l2->p4.Pt() > 15.))
    return {};

  return std::make_tuple(leptonCat, l1, l2);
}

void ZGammaTrees::FillMoreVariables(std::vector<Jet> const &jets) {
  run_ = *srcRun_;
  lumi_ = *srcLumi_;
  event_ = *srcEvent_;
  jetSize_ = std::min<int>(jets.size(), maxSize_);
}
