#include <TrileptonTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TLorentzVector.h>
#include <TVector2.h>

#include <Utils.h>


namespace po = boost::program_options;


int const TrileptonTrees::maxSize_;


TrileptonTrees::TrileptonTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      ptMissCut_{options.GetAs<double>("ptmiss-cut")},
      triggerFilter_{dataset, options, &runSampler_},
      srcRun_{dataset.Reader(), "run"},
      srcLumi_{dataset.Reader(), "luminosityBlock"},
      srcEvent_{dataset.Reader(), "event"},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  if (isSim_) {
    srcLHEVpt_.reset(new TTreeReaderValue<Float_t>(dataset.Reader(), "LHE_Vpt"));

    auto const &node = dataset.Info().Parameters()["zz_2l2nu"];

    if (node and not node.IsNull() and node.as<bool>())
      genZZBuilder_.emplace(dataset);
  }

  tauBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});

  CreateWeightBranches();

  AddBranch("lepton_cat", &leptonCat_);
  AddBranch("jet_cat", &jetCat_);
  AddBranch("ll_pt", &llPt_);
  AddBranch("ll_eta", &llEta_);
  AddBranch("ll_phi", &llPhi_);
  AddBranch("ll_mass", &llMass_);
  AddBranch("l3_pt", &l3Pt_);
  AddBranch("l3_eta", &l3Eta_);
  AddBranch("l3_phi", &l3Phi_);
  AddBranch("l3_mass", &l3Mass_);
  AddBranch("mTw", &mTw_);
  AddBranch("ptmiss", &missPt_);
  AddBranch("ptmiss_phi", &missPhi_);
  AddBranch("ptmiss_significance", &missSignificance_);
  AddBranch("ptmiss_significance_corrected", &missSignificanceCorrected_);
  AddBranch("dphi_visibles_ptmiss", &dPhiVisiblesPtmiss_);
  AddBranch("mT", &mT_);
  AddBranch("num_pv_good", &numPVGood_);
  AddBranch("dijet_mass", &dijetMass_);

  if (storeMoreVariables_) {
    AddBranch("run", &run_);
    AddBranch("lumi", &lumi_);
    AddBranch("event", &event_);

    if (genZZBuilder_)
      AddBranch("gen_mzz", &genMZZ_);

    AddBranch("lepton_charge", leptonCharge_, "lepton_charge[3]/I");
    AddBranch("lepton_pt", leptonPt_, "lepton_pt[3]/F");
    AddBranch("lepton_eta", leptonEta_, "lepton_eta[3]/F");
    AddBranch("lepton_phi", leptonPhi_, "lepton_phi[3]/F");
    AddBranch("lepton_mass", leptonMass_, "lepton_mass[3]/F");

    AddBranch("jet_size", &jetSize_);
    AddBranch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    AddBranch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    AddBranch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    AddBranch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }

  datasetLHEVptUpperLimitInc_.reset();
  auto const LHEVptUpperLimitIncSettingsNode = dataset.Info().Parameters()["LHE_Vpt_upper_limit_inc"];
  if (LHEVptUpperLimitIncSettingsNode and not LHEVptUpperLimitIncSettingsNode.IsNull()) {
    datasetLHEVptUpperLimitInc_.emplace(LHEVptUpperLimitIncSettingsNode.as<Float_t>());
  }
}


po::options_description TrileptonTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  optionsDescription.add_options()
    ("ptmiss-cut", po::value<double>()->default_value(30.),
     "Minimal missing pt");
  return optionsDescription;
}


bool TrileptonTrees::ProcessEvent() {
  if (datasetLHEVptUpperLimitInc_.has_value() and not (*srcLHEVpt_->Get() <= datasetLHEVptUpperLimitInc_.value()))
    return false;

  if (not ApplyCommonFilters())
    return false;

  auto const leptonResult = CheckLeptons();
  if (not leptonResult)
    return false;

  if (tauBuilder_.Get().size() > 0)
    return false;

  auto const &[leptonCat, l1, l2, l3] = leptonResult.value();
  // switch (leptonCat) {
  //   case LeptonCat::kEE:
  //     if (not triggerFilter_("ee"))
  //       return false;
  //     break;
  //   case LeptonCat::kMuMu:
  //     if (not triggerFilter_("mumu"))
  //       return false;
  //     break;
  //   case LeptonCat::kEMu:
  //     if (not triggerFilter_("emu"))
  //       return false;
  //     break;
  // }

  leptonCat_ = int(leptonCat);
  TLorentzVector const p4LL = l1->p4 + l2->p4;
  llPt_ = p4LL.Pt();
  llEta_ = p4LL.Eta();
  llPhi_ = p4LL.Phi();
  llMass_ = p4LL.M();

  l3Pt_ = l3->p4.Pt();
  l3Eta_ = l3->p4.Eta();
  l3Phi_ = l3->p4.Phi();
  l3Mass_ = l3->p4.M();

  if (not (std::abs(p4LL.M() - kNominalMZ_) < zMassWindow_))
    return false;

//  if (not (p4LL.Pt() > minPtLL_))
//    return false;


  auto const &p4Miss = ptMissBuilder_.Get().p4;
  missPt_ = p4Miss.Pt();
  missPhi_ = p4Miss.Phi();
  missSignificance_ = ptMissBuilder_.Get().significance;
  missSignificanceCorrected_ = ptMissBuilder_.Get().significance_corrected;

  mTw_ = std::sqrt(2*l3Pt_*missPt_*(1-std::cos(l3Phi_ - missPhi_)));
//  if (mTw_ < mTwCut_) //Cut value set in .h file
//    return false;

  if (p4Miss.Pt() < ptMissCut_)
    return false;

     // ptmiss_significance_corrected > 3 cut
  if (missSignificanceCorrected_ <= 3)
    return false;

//  if (std::abs(
//        TVector2::Phi_mpi_pi(p4LL.Phi() - p4Miss.Phi())) < minDphiLLPtMiss_)
//    return false;


  auto const &jets = jetBuilder_.Get();

  for (auto const &jet : jets) {
    if (bTagger_(jet))
      return false;

//    if (std::abs(TVector2::Phi_mpi_pi(
//            jet.p4.Phi() - p4Miss.Phi())) < minDphiJetPtMiss_)
//      return false;
  }

  dPhiVisiblesPtmiss_ = DPhiPtMiss({&jetBuilder_, &muonBuilder_, &electronBuilder_});

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (jets.size() == 1)
    jetCat_ = int(JetCat::kEq1J);
  else
    jetCat_ = int(JetCat::kGEq2J);


  if (jets.size() >= 2) {
    auto dijet_p4 =  jets[0].p4 + jets[1].p4;
    dijetMass_ = dijet_p4.M();
  } else {
    dijetMass_ = -1.;
  }

  double const eT =
      std::sqrt(std::pow(p4LL.Pt(), 2) + std::pow(p4LL.M(), 2))
      + std::sqrt(std::pow(p4Miss.Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(std::pow(eT, 2) - std::pow((p4LL + p4Miss).Pt(), 2));

  numPVGood_ = *srcNumPVGood_;


  if (storeMoreVariables_)
    FillMoreVariables({*l1, *l2, *l3}, jets);

  FillTree();
  return true;
}


std::optional<std::tuple<TrileptonTrees::LeptonCat, Lepton const *, Lepton const *, Lepton const *>>
TrileptonTrees::CheckLeptons() const {
  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  if (looseElectrons.size() + looseMuons.size() != 3)
    return {};

  LeptonCat leptonCat;
  Lepton const *l1 = nullptr, *l2 = nullptr, *l3 = nullptr;

  if (tightElectrons.size() == 3) {
    leptonCat = LeptonCat::kEEE;
    Float_t M_OSSF_residual = 9999;
    for (unsigned ei=0; ei<3; ei++){
      for (unsigned ej=ei+1; ej<3; ej++){
        Lepton const *tmpe1, *tmpe2;
        TLorentzVector tmpZee;
        tmpe1 = &tightElectrons[ei];
        tmpe2 = &tightElectrons[ej];
        if (tmpe1->charge * tmpe2->charge > 0) continue;
        tmpZee = tmpe1->p4 + tmpe2->p4;
        if(abs(tmpZee.M()-kNominalMZ_) < M_OSSF_residual){
          M_OSSF_residual = abs(tmpZee.M()-kNominalMZ_);
          l1 = &tightElectrons[ei];
          l2 = &tightElectrons[ej];
          l3 = &tightElectrons[3-ei-ej];
        }
      }
    }
    if(M_OSSF_residual == 9999) return {};
  }
  else if(tightMuons.size() == 3){
    leptonCat = LeptonCat::kMuMuMu;
    Float_t M_OSSF_residual = 9999;
    for (unsigned mi=0; mi<3; mi++){
      for (unsigned mj=mi+1; mj<3; mj++){
        Lepton const *tmpm1, *tmpm2;
        TLorentzVector tmpZmm;
        tmpm1 = &tightMuons[mi];
        tmpm2 = &tightMuons[mj];
        if (tmpm1->charge * tmpm2->charge > 0) continue;
        tmpZmm = tmpm1->p4 + tmpm2->p4;
        if(abs(tmpZmm.M()-kNominalMZ_) < M_OSSF_residual){
          M_OSSF_residual = abs(tmpZmm.M()-kNominalMZ_);
          l1 = &tightMuons[mi];
          l2 = &tightMuons[mj];
          l3 = &tightMuons[3-mi-mj];
        }
      }
    }
    if(M_OSSF_residual == 9999) return {};
  }
  else if ((tightElectrons.size() == 2) and (tightMuons.size() == 1)) {
    leptonCat = LeptonCat::kEEMu;
    l1 = &tightElectrons[0];
    l2 = &tightElectrons[1];
    l3 = &tightMuons[0];
  }
  else if ((tightElectrons.size() == 1) and (tightMuons.size() == 2)) {
    leptonCat = LeptonCat::kMUMUE;
    l1 = &tightMuons[0];
    l2 = &tightMuons[1];
    l3 = &tightElectrons[0];
  }else
    return {};

  if (not (l1->charge * l2->charge < 0))
    return {};


  if (l1->p4.Pt() < l2->p4.Pt())
    std::swap(l1, l2);

  if (not (l1->p4.Pt() > 30. && l2->p4.Pt() > 20. && l3->p4.Pt() > 20. ))
    return {};

  return std::make_tuple(leptonCat, l1, l2, l3);
}


void TrileptonTrees::FillMoreVariables(
    std::array<Lepton, 3> const &leptons, std::vector<Jet> const &jets) {

  run_ = *srcRun_;
  lumi_ = *srcLumi_;
  event_ = *srcEvent_;

  if (genZZBuilder_)
    genMZZ_ = genZZBuilder_->P4ZZ().M();

  for (int i = 0; i < 3; ++i) {
    leptonCharge_[i] = leptons[i].charge;
    auto const &p4 = leptons[i].p4;
    leptonPt_[i] = p4.Pt();
    leptonEta_[i] = p4.Eta();
    leptonPhi_[i] = p4.Phi();
    leptonMass_[i] = p4.M();
  }

  jetSize_ = std::min<int>(jets.size(), maxSize_);

  for (int i = 0; i < jetSize_; ++i) {
    auto const &p4 = jets[i].p4;
    jetPt_[i] = p4.Pt();
    jetEta_[i] = p4.Eta();
    jetPhi_[i] = p4.Phi();
    jetMass_[i] = p4.M();
  }
}
