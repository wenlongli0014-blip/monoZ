#include <JetBuilder.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include <Utils.h>


JetBuilder::JetBuilder(
    Dataset &dataset, Options const &options, TabulatedRngEngine &rngEngine,
    PileUpIdFilter const *pileUpIdFilter)
    : CollectionBuilder{dataset.Reader()},
      genJetBuilder_{nullptr}, pileUpIdFilter_{pileUpIdFilter},
      minPtType1Corr_{15.}, ptMissEeNoise_{false}, ptMissPogJets_{false},
      isSim_{dataset.Info().IsSimulation()},
      jetCorrector_{dataset, options, rngEngine},
      srcPt_{dataset.Reader(), "Jet_pt"},
      srcEta_{dataset.Reader(), "Jet_eta"},
      srcPhi_{dataset.Reader(), "Jet_phi"},
      srcMass_{dataset.Reader(), "Jet_mass"},
      srcArea_{dataset.Reader(), "Jet_area"},
      srcRawFactor_{dataset.Reader(), "Jet_rawFactor"},
      srcChEmEF_{dataset.Reader(), "Jet_chEmEF"},
      srcNeEmEF_{dataset.Reader(), "Jet_neEmEF"},
      srcMuonFraction_{dataset.Reader(), "Jet_muonSubtrFactor"},
      srcBTag_{dataset.Reader(), (Options::NodeAs<std::string>(
        options.GetConfig(), {"b_tagger", "branch_name"})).c_str()},
      srcId_{dataset.Reader(), "Jet_jetId"},
      srcPileUpId_{dataset.Reader(), "Jet_puId"},
      puRho_{dataset.Reader(), "fixedGridRhoFastjetAll"},
      softRawPt_{dataset.Reader(), "CorrT1METJet_rawPt"},
      softEta_{dataset.Reader(), "CorrT1METJet_eta"},
      softPhi_{dataset.Reader(), "CorrT1METJet_phi"},
      softArea_{dataset.Reader(), "CorrT1METJet_area"},
      softMuonFraction_{dataset.Reader(), "CorrT1METJet_muonSubtrFactor"} {

  auto const jetConfig = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"jets"});
  jetIdBit_ = Options::NodeAs<int>(jetConfig, {"jet_id_bit"});
  minPt_ = Options::NodeAs<double>(jetConfig, {"min_pt"});
  maxAbsEta_ = Options::NodeAs<double>(jetConfig, {"max_abs_eta"});

  if (pileUpIdFilter_)
    LOG_DEBUG << "PileUpIdFilter is registered in JetBuilder.";
  if (pileUpIdFilter_) {
    std::tie(pileUpIdMinPt_, pileUpIdMaxPt_) = pileUpIdFilter_->GetPtRange();
  } else {
    // Defaults from https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID?rev=61#Recommendations_for_13_TeV_data
    pileUpIdMinPt_ = 15.;
    pileUpIdMaxPt_ = 50.;
  }

  auto const ptMissConfig = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"ptmiss"});
  ptMissJer_ = Options::NodeAs<bool>(ptMissConfig, {"jer"});
  if (auto const node = ptMissConfig["pog_jets"]; node)
    ptMissPogJets_ = node.as<bool>();
  if (auto const node = ptMissConfig["ptmiss_fix_ee_2017"]; node)
    ptMissEeNoise_ = node.as<bool>();

  if (isSim_) {
    srcHadronFlavour_.emplace(dataset.Reader(), "Jet_hadronFlavour");
    srcPartonFlavour_.emplace(dataset.Reader(), "Jet_partonFlavour");
    srcGenJetIdx_.emplace(dataset.Reader(), "Jet_genJetIdx");
  }
}


std::vector<Jet> const &JetBuilder::Get() const {
  Update();
  return jets_;
}


std::vector<Jet> const &JetBuilder::GetRejected() const {
  Update();
  return rejectedJets_;
}

std::vector<Jet> const &JetBuilder::GetLowPt() const {
  Update();
  return lowptJets_;
}

void JetBuilder::SetGenJetBuilder(GenJetBuilder const *genJetBuilder) {
  if (isSim_)
    genJetBuilder_ = genJetBuilder;
}


void JetBuilder::AddType1Correction(
    TLorentzVector const &rawP4, double area,
    double jecOrig, double jecNew, double jerFactor,
    double emFraction, double muonFraction) const {
  // If jets are to be treated as in the standard type 1 correction [1], skip
  // those with a high EM energy fraction and subtract muons from jets. The
  // muon subtraction is approximate in that it's implemented as a rescaling of
  // the momentum. When JEC are applied below, they effectively rescale also the
  // muon's momenta; this is how it's done in the standard procedure.
  // [1] https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/issues/49
  if (ptMissPogJets_ and emFraction > 0.9)
    return;
  double const rescale = (ptMissPogJets_) ? 1. - muonFraction : 1.;

  double corrFactorNew = jecNew;
  if (ptMissJer_)
    corrFactorNew *= jerFactor;

  // Check if this jet should be subjected to the EE noise mitigation procedure
  bool isEeNoise = false;
  if (ptMissEeNoise_) {
    double const absEta = std::abs(rawP4.Eta());
    // Thresholds are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription?rev=91#Instructions_for_2017_data_with
    isEeNoise = (rawP4.Pt() < 50. and absEta > 2.65 and absEta < 3.139);
  }

  // The normal type 1 correction for jets not affected by the EE noise
  double const jecL1 = jetCorrector_.GetJecL1(rawP4, area);
  if (not isEeNoise and rawP4.Pt() * corrFactorNew > minPtType1Corr_)
    AddMomentumShift(rawP4 * jecL1 * rescale, rawP4 * corrFactorNew * rescale);

  // If the EE noise mitigation is enabled, the computation starts from an
  // adjusted ptmiss instead of the raw one, and some of the contributions in it
  // have to be removed. Technically, the removal procedure is equivalent to
  // applying the type 1 correction with affected jets but using the same JEC as
  // when NanoAOD was produced. See [1] for details.
  // [1] https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/issues/49
  if (isEeNoise) {
    // Always apply the filtering and rescaling as if ptMissPogJets_ is true in
    // order to match closer the compuatation of the contribution that is being
    // removed here
    if (emFraction > 0.9)
      return;
    double const rescale = 1. - muonFraction;

    if (rawP4.Pt() * jecOrig > minPtType1Corr_)
      AddMomentumShift(rawP4 * jecL1 * rescale, rawP4 * jecOrig * rescale);
  }
}


void JetBuilder::Build() const {
  jetCorrector_.UpdateIov();
  ProcessJets();

  // Soft jets not included into the main collection contribute to the type 1
  // correction of missing pt nonetheless. Account for them.
  ProcessSoftJets();
}


GenJet const *JetBuilder::FindGenMatch(TLorentzVector const &p4,
                                       double ptResolution) const {
  if (not genJetBuilder_)
    return nullptr;

  // Find the closest generator-level jet in the (eta, phi) metric, with an
  // additional requirement that the difference in pt is loosely compatible with
  // the resolution. The angular distance must not exceed a half of jet radius.
  // These requirements are taken from [1].
  // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
  GenJet const *match = nullptr;
  double const maxDPt = ptResolution * p4.Pt() * 3;  // Use 3 sigma
  double curMinDR2 = std::pow(0.4 / 2, 2);  // Use half of jet radius

  for (auto const &genJet : genJetBuilder_->Get()) {
    double const dR2 = utils::DeltaR2(p4, genJet.p4);

    if (dR2 < curMinDR2 and std::abs(p4.Pt() - genJet.p4.Pt()) < maxDPt) {
      match = &genJet;
      curMinDR2 = dR2;
    }
  }

  return match;
}


double JetBuilder::GetJerFactor(
    TLorentzVector const &corrP4, int rngChannel) const {
  double const ptResolution = jetCorrector_.GetPtResolution(corrP4);
  GenJet const *genJet = FindGenMatch(corrP4, ptResolution);
  return jetCorrector_.GetJerFactor(
      corrP4, genJet, ptResolution, rngChannel);
}


void JetBuilder::ProcessJets() const {
  jets_.clear();
  lowptJets_.clear();
  rejectedJets_.clear();

 for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    Jet jet;
    jet.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    // double const jecRaw = 1. / (1 - srcRawFactor_[i]);
    // TLorentzVector const rawP4 = jet.p4  * (1. / jecRaw);
    // double const area = srcArea_[i];
    // double const jecNominal = jetCorrector_.GetJecFull(rawP4, area);
    double const jecNominal = 1. / (1 - srcRawFactor_[i]);
    TLorentzVector const rawP4 = jet.p4  * (1. / jecNominal);
    double jecUncFactor = 1., jerFactor = 1.;
    if (isSim_) {
      jecUncFactor = jetCorrector_.GetJecUncFactor(jet.p4);
      jerFactor = GetJerFactor(jet.p4, i);
    }
    jet.p4 *= jecUncFactor * jerFactor;

    AddType1Correction(
        rawP4, srcArea_[i], jecNominal, jecNominal * jecUncFactor, jerFactor,
        srcChEmEF_[i] + srcNeEmEF_[i], srcMuonFraction_[i]);

    // Kinematical cuts and ID selection for jets to be stored in the collection
    if (not (srcId_[i] & 1 << jetIdBit_))
      continue;
    if (std::abs(jet.p4.Eta()) > maxAbsEta_)
      continue;

    // Perform angular cleaning
    if (IsDuplicate(jet.p4, 0.4))
      continue;

    jet.bTag = srcBTag_[i];
    if (isSim_)
      jet.SetFlavours(srcHadronFlavour_->At(i), srcPartonFlavour_->At(i));

    bool const puIdAccepted = SetPileUpInfo(jet, i);
    if(jet.p4.Pt() > minPt_)
    {	
      if (puIdAccepted)
        jets_.emplace_back(jet);
      else
        rejectedJets_.emplace_back(jet);
    }
    else if (jet.p4.Pt() > 20.0 && jet.p4.Pt() < minPt_)
      lowptJets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
  std::sort(rejectedJets_.begin(), rejectedJets_.end(), PtOrdered);
  std::sort(lowptJets_.begin(), lowptJets_.end(), PtOrdered);
}


void JetBuilder::ProcessSoftJets() const {
  for (int i = 0; i < int(softRawPt_.GetSize()); ++i) {
    // Jet energy is not stored, but it's not used for missing pt. Set the mass
    // to 0.
    TLorentzVector rawP4;
    rawP4.SetPtEtaPhiM(softRawPt_[i], softEta_[i], softPhi_[i], 0.);

    double const area = softArea_[i];
    double const jecNominal = jetCorrector_.GetJecFull(rawP4, area);
    double jecFull = jecNominal, jerFactor = 1.;

    if (isSim_) {
      TLorentzVector const corrP4{rawP4 * jecNominal};
      jecFull *= jetCorrector_.GetJecUncFactor(corrP4);
      if (ptMissJer_)
        jerFactor = GetJerFactor(corrP4, srcPt_.GetSize() + i);
    }

    AddType1Correction(
        rawP4, area, jecNominal, jecFull, jerFactor,
        0. /* EM fractions are not stored */, softMuonFraction_[i]);
  }
}


bool JetBuilder::SetPileUpInfo(Jet &jet, int index) const {
  if (jet.p4.Pt() < pileUpIdMinPt_ or jet.p4.Pt() > pileUpIdMaxPt_) {
    jet.pileUpId = Jet::PileUpId::PassThrough;
  } else {
    int const id = srcPileUpId_[index];
    if (id & 1)
      jet.pileUpId = Jet::PileUpId::Tight;
    else if (id & 1 << 1)
      jet.pileUpId = Jet::PileUpId::Medium;
    else if (id & 1 << 2)
      jet.pileUpId = Jet::PileUpId::Loose;
    else
      jet.pileUpId = Jet::PileUpId::None;
  }

  if (isSim_) {
    // The standard matching [1] is done based on dR, with the cut-off dR < 0.4.
    // However, each particle-level jet can be matched to one reconstructed jet
    // at maximum [2] to guarantee an unambiguous matching.
    // [1] https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_22/PhysicsTools/PatAlgos/python/mcMatchLayer0/jetMatch_cfi.py#L18-L28
    // [2] https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_22/CommonTools/UtilAlgos/interface/PhysObjectMatcher.h#L150-L159
    jet.isPileUp = (srcGenJetIdx_->At(index) == -1);
  }

  if (pileUpIdFilter_ and not (*pileUpIdFilter_)(jet)) {
    LOG_TRACE << "Pileup ID filter rejects jets with pt " << jet.p4.Pt()
        << " GeV, eta " << jet.p4.Eta() << ", and pileup ID WP "
        << int(jet.pileUpId) << ".";
    return false;
  }

  return true;
}
