#include <NrbTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TLorentzVector.h>
#include <TVector2.h>

#include <Utils.h>


namespace po = boost::program_options;


int const NrbTrees::maxSize_;


NrbTrees::NrbTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      ptMissCut_{options.GetAs<double>("ptmiss-cut")},
      triggerFilter_{dataset, options, &runSampler_},
      leptonEff_{dataset, options, &electronBuilder_, &muonBuilder_, isSim_? 2 : 1},
      triggerEff_{dataset, options, &electronBuilder_, &muonBuilder_, isSim_? 2 : 1},
      srcEvent_{dataset.Reader(), "event"},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  if (isSim_) {
    auto const &node = dataset.Info().Parameters()["zz_2l2nu"];

    if (node and not node.IsNull() and node.as<bool>())
      genZZBuilder_.emplace(dataset);
  }

  CreateWeightBranches();

  AddBranch("lepton_cat", &leptonCat_);
  AddBranch("jet_cat", &jetCat_);
  AddBranch("ll_pt", &llPt_);
  AddBranch("ll_eta", &llEta_);
  AddBranch("ll_phi", &llPhi_);
  AddBranch("ll_mass", &llMass_);
  AddBranch("ptmiss", &missPt_);
  AddBranch("ptmiss_phi", &missPhi_);
  AddBranch("mT", &mT_);
  AddBranch("num_pv_good", &numPVGood_);
  AddBranch("l1_pt", &l1Pt_);
  AddBranch("l2_pt", &l1Pt_);
  AddBranch("l1_eta", &l2Eta_);
  AddBranch("l2_eta", &l2Eta_);
  AddBranch("btag_loose", &btagLoose_);
  AddBranch("btag_medium", &btagMedium_);
  AddBranch("btag_tight", &btagTight_);
  AddBranch("btag_loose_lowpt", &btagLooseLowPt_);
  AddBranch("btag_medium_lowpt", &btagMediumLowPt_);
  AddBranch("btag_tight_lowpt", &btagTightLowPt_);

  AddBranch("tmp_weight", &weight_);
  int const triggerVariations = triggerEff_.NumVariations();
  triggerWeights_.resize(triggerVariations);
  for (int i = 0; i < triggerVariations; ++i) {
    auto const name = "tmp_weight_"
        + std::string{triggerEff_.VariationName(i)};
    AddBranch(name.c_str(), &triggerWeights_[i]);
  }

  int const leptonVariations = 8;//leptonEff_.NumVariations();
  leptonWeights_.resize(leptonVariations);
  for (int i = 0; i < leptonVariations; ++i) {
    auto const name = "tmp_weight_"
        + std::string{leptonEff_.VariationName(i < 4 ? i : i + 4)};
    AddBranch(name.c_str(), &leptonWeights_[i]);
  }

  if (storeMoreVariables_) {
    AddBranch("event", &event_);

    if (genZZBuilder_)
      AddBranch("gen_mzz", &genMZZ_);

    AddBranch("lepton_charge", leptonCharge_, "lepton_charge[2]/I");
    AddBranch("lepton_pt", leptonPt_, "lepton_pt[2]/F");
    AddBranch("lepton_eta", leptonEta_, "lepton_eta[2]/F");
    AddBranch("lepton_phi", leptonPhi_, "lepton_phi[2]/F");
    AddBranch("lepton_mass", leptonMass_, "lepton_mass[2]/F");

    AddBranch("jet_size", &jetSize_);
    AddBranch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    AddBranch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    AddBranch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    AddBranch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }
}


po::options_description NrbTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  optionsDescription.add_options()
    ("ptmiss-cut", po::value<double>()->default_value(50.), 
     "Minimal missing pt");
  return optionsDescription;
}


bool NrbTrees::ProcessEvent() {
  if (not ApplyCommonFilters())
    return false;

  auto const leptonResult = CheckLeptons();
  if (not leptonResult)
    return false;

  if (tauBuilder_.Get().size() > 0)
    return false;

  auto const &[leptonCat, l1, l2] = leptonResult.value();
  switch (leptonCat) {
    case LeptonCat::kEE:
      if (not triggerFilter_("ee"))
        return false;
      break;
    case LeptonCat::kMuMu:
      if (not triggerFilter_("mumu"))
        return false;
      break;
    case LeptonCat::kEMu:
      if (not triggerFilter_("emu"))
        return false;
      break;
  }

  leptonCat_ = int(leptonCat);
  l1Pt_ = l1->p4.Pt();
  l2Pt_ = l2->p4.Pt();
  l1Eta_ = l1->p4.Eta();
  l2Eta_ = l1->p4.Eta();
  TLorentzVector const p4LL = l1->p4 + l2->p4;
  llPt_ = p4LL.Pt();
  llEta_ = p4LL.Eta();
  llPhi_ = p4LL.Phi();
  llMass_ = p4LL.M();

  if (p4LL.M() < 50.)
    return false;

  if (p4LL.Pt() < 25.)
    return false;


  auto const &p4Miss = ptMissBuilder_.Get().p4;
  missPt_ = p4Miss.Pt();
  missPhi_ = p4Miss.Phi();

  if (p4Miss.Pt() < ptMissCut_)
    return false;

  if (std::abs(
        TVector2::Phi_mpi_pi(p4LL.Phi() - p4Miss.Phi())) < minDphiLLPtMiss_)
    return false;


  auto const &jets = jetBuilder_.Get();

  btagLoose_ = false;
  btagMedium_ = false;
  btagTight_ = false;
  
  for (auto const &jet : jets) {
    if (std::abs(TVector2::Phi_mpi_pi(
            jet.p4.Phi() - p4Miss.Phi())) < minDphiJetPtMiss_)
      return false;
    if (bTagger_.Loose(jet)) btagLoose_ = true;
    if (bTagger_.Medium(jet)) btagMedium_ = true;
    if (bTagger_.Tight(jet)) btagTight_ = true;
  }

  auto const &jetsL = jetBuilder_.GetLowPt();

  btagLooseLowPt_ = false;
  btagMediumLowPt_ = false;
  btagTightLowPt_ = false;

  for (auto const &jetl : jetsL) {
    if (bTagger_.Loose(jetl)) btagLooseLowPt_ = true;
    if (bTagger_.Medium(jetl)) btagMediumLowPt_ = true;
    if (bTagger_.Tight(jetl)) btagTightLowPt_ = true;
  }
  // if (DPhiPtMiss({&jetBuilder_, &muonBuilder_, &electronBuilder_})
  //     < minDphiLeptonsJetsPtMiss_)
  //   return false;

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (jets.size() == 1)
    jetCat_ = int(JetCat::kEq1J);
  else
    jetCat_ = int(JetCat::kGEq2J);


  double const eT =
      std::sqrt(std::pow(p4LL.Pt(), 2) + std::pow(p4LL.M(), 2))
      + std::sqrt(std::pow(p4Miss.Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(std::pow(eT, 2) - std::pow((p4LL + p4Miss).Pt(), 2));

  numPVGood_ = *srcNumPVGood_;


  // for nominal weight
  if (leptonCat == LeptonCat::kEMu) 
    weight_ = leptonEff_.TransferFactor(l1, l2, 0) 
      * triggerEff_.TransferFactor(l1, l2, 0);
  else 
    weight_ = 0.; 
  // for systematics: muonEff_syst_up&down, muonEff_stat_up&down
  // electronEff_syst_up&down, electronEff_stat_up&down
  for (int i = 0; i < 8; i++) {
    if (leptonCat == LeptonCat::kEMu) 
      leptonWeights_[i] = leptonEff_.TransferFactor(l1, l2, i < 4 ? i + 1 : i + 5) *
        triggerEff_.TransferFactor(l1, l2, 0);
    else 
      leptonWeights_[i] = 0.;
  }	
  // for systematics: triggerEE_up&down, triggerMuMu_up&down, triggerEMu_up&down
  for (int i = 0; i < triggerEff_.NumVariations(); i++) {
    if (leptonCat == LeptonCat::kEMu) 
      triggerWeights_[i] = leptonEff_.TransferFactor(l1, l2, 0)
        *triggerEff_.TransferFactor(l1, l2, i + 1);
    else 
      triggerWeights_[i] = 0.;
  }	
  if (storeMoreVariables_)
    FillMoreVariables({*l1, *l2}, jets);

  FillTree();
  return true;
}


std::optional<std::tuple<NrbTrees::LeptonCat, Lepton const *, Lepton const *>>
NrbTrees::CheckLeptons() const {
  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  if (looseElectrons.size() + looseMuons.size() != 2)
    return {};

  LeptonCat leptonCat;
  Lepton const *l1, *l2;

  if (tightElectrons.size() == 2) {
    leptonCat = LeptonCat::kEE;
    l1 = &tightElectrons[0];
    l2 = &tightElectrons[1];
    if (l1->p4.Pt() < l2->p4.Pt())
      std::swap(l1, l2);
  } else if (tightMuons.size() == 2) {
    leptonCat = LeptonCat::kMuMu;
    l1 = &tightMuons[0];
    l2 = &tightMuons[1];
    if (l1->p4.Pt() < l2->p4.Pt())
      std::swap(l1, l2);
  } else if (tightElectrons.size() == 1 and tightMuons.size() == 1) {
    leptonCat = LeptonCat::kEMu;
    l1 = &tightElectrons[0];
    l2 = &tightMuons[0];
  } else
    return {};

  return std::make_tuple(leptonCat, l1, l2);
}


void NrbTrees::FillMoreVariables(
    std::array<Lepton, 2> const &leptons, std::vector<Jet> const &jets) {

  event_ = *srcEvent_;

  if (genZZBuilder_)
    genMZZ_ = genZZBuilder_->P4ZZ().M();

  for (int i = 0; i < 2; ++i) {
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
