#include <PtMissBuilder.h>

#include <Logger.h>
#include <MetXYCorrections.h>


PtMissBuilder::PtMissBuilder(Dataset &dataset, Options const &options)
    : syst_{Syst::None}, applyEeNoiseMitigation_{false},
      cache_{dataset.Reader()},
      isSim_{dataset.Info().IsSimulation()},
      srcNumPV_{dataset.Reader(), "PV_npvs"},
      srcPt_{dataset.Reader(), "RawMET_pt"},
      srcPhi_{dataset.Reader(), "RawMET_phi"},
      srcRun_{dataset.Reader(), "run"},
      srcCOVxx_{dataset.Reader(), "MET_covXX"},
      srcCOVxy_{dataset.Reader(), "MET_covXY"},
      srcCOVyy_{dataset.Reader(), "MET_covYY"} {

  auto const config = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"ptmiss"});
  if (auto const node = config["fix_ee_2017"]; node and node.as<bool>()) {
    applyEeNoiseMitigation_ = true;
    LOG_DEBUG << "Will apply EE noise mitigation in missing pt.";
  }
  metXYCorrectionYear_ = config["XY_corrections"].as<std::string>();

  std::vector<std::string> listOfGoodYears;

  if (auto const node = config["is_ul"]; node and node.as<bool>()) {
    isUL_ = true;
    if (!isSim_ && (metXYCorrectionYear_ == "2016APV" || metXYCorrectionYear_ == "2016nonAPV")) {
      metXYCorrectionYear_ = "2016";
    }
    listOfGoodYears = {"2016", "2016APV", "2016nonAPV", "2017", "2018"};
    LOG_DEBUG << "Will choose MET XY corrections for UL datasets year " << metXYCorrectionYear_;
  } else {
    isUL_ = false;
    listOfGoodYears = {"2016", "2017", "2018"};
    LOG_DEBUG << "Will choose MET XY corrections for pre-UL datasets year " << metXYCorrectionYear_;
  }

  applyXYCorrections_ = (std::find(std::begin(listOfGoodYears),
    std::end(listOfGoodYears), metXYCorrectionYear_) !=
    std::end(listOfGoodYears));

  if (!applyXYCorrections_)
      LOG_DEBUG << "MET XY corrections will NOT be applied.";

  if (applyEeNoiseMitigation_) {
    srcDefaultPt_.emplace(dataset.Reader(), "MET_pt");
    srcDefaultPhi_.emplace(dataset.Reader(), "MET_phi");
    srcFixedPt_.emplace(dataset.Reader(), "METFixEE2017_pt");
    srcFixedPhi_.emplace(dataset.Reader(), "METFixEE2017_phi");
    srcSignificance_.emplace(dataset.Reader(), "METFixEE2017_significance");
  } else {
    srcSignificance_.emplace(dataset.Reader(), "MET_significance");
  }

  std::string const systLabel{options.GetAs<std::string>("syst")};
  if (systLabel == "metuncl_up")
    syst_ = Syst::UnclEnergyUp;
  else if (systLabel == "metuncl_down")
    syst_ = Syst::UnclEnergyDown;
  else
    syst_ = Syst::None;

  if (syst_ != Syst::None) {
    std::string const name{(applyEeNoiseMitigation_) ? "METFixEE2017" : "MET"};
    srcUnclEnergyUpDeltaX_.emplace(
        dataset.Reader(), (name + "_MetUnclustEnUpDeltaX").c_str());
    srcUnclEnergyUpDeltaY_.emplace(
        dataset.Reader(), (name + "_MetUnclustEnUpDeltaY").c_str());
    LOG_DEBUG << "Will apply a variation in unclustered momentum in ptmiss.";
  }
}


PtMiss const &PtMissBuilder::Get() const {
  if (cache_.IsUpdated())
    Build();

  return ptMiss_;
}


void PtMissBuilder::PullCalibration(
    std::initializer_list<CollectionBuilderBase const *> builders) {
  for (auto const *b : builders)
    calibratingBuilders_.emplace_back(b);
}


void PtMissBuilder::Build() const {
  ptMiss_.p4.SetPtEtaPhiM(*srcPt_, 0., *srcPhi_, 0.);
  ptMiss_.significance = **srcSignificance_;
  //std::cout<<"raw ptmiss pt:"<<*srcPt_ << "  phi:"<<*srcPhi_<<std::endl;
  if (applyEeNoiseMitigation_) {
    // Add (METFixEE2017 - MET) to the raw missing pt. This difference PF
    // candidates in the EE, both unclustered and also included in soft jets.
    // It also picks up a part of the type 1 correction from MET, which is
    // removed in JetBuilder. See [1] for details.
    // [1] https://hypernews.cern.ch/HyperNews/CMS/get/met/710.html
    TLorentzVector p4;
    p4.SetPtEtaPhiM(**srcFixedPt_, 0., **srcFixedPhi_, 0.);
    ptMiss_.p4 += p4;
    p4.SetPtEtaPhiM(**srcDefaultPt_, 0., **srcDefaultPhi_, 0.);
    ptMiss_.p4 -= p4;
  }

  for (auto const *builder : calibratingBuilders_)
    ptMiss_.p4 -= builder->GetSumMomentumShift();

  // Apply MET XY corrections, if the corresponding option has been set.
  if (applyXYCorrections_){
    std::pair<double, double> corrected_met_metPhi = 
        metXYCorrections::METXYCorr_Met_MetPhi(ptMiss_.p4.Pt(),
        ptMiss_.p4.Phi(), *srcRun_, metXYCorrectionYear_, isSim_,
        *srcNumPV_, isUL_, false);
    ptMiss_.p4.SetPtEtaPhiM(
        corrected_met_metPhi.first, 0, corrected_met_metPhi.second, 0);
  }

  if (syst_ == Syst::UnclEnergyUp) {
    ptMiss_.p4.SetPx(ptMiss_.p4.Px() + **srcUnclEnergyUpDeltaX_);
    ptMiss_.p4.SetPy(ptMiss_.p4.Py() + **srcUnclEnergyUpDeltaY_);
  } else if (syst_ == Syst::UnclEnergyDown) {
    ptMiss_.p4.SetPx(ptMiss_.p4.Px() - **srcUnclEnergyUpDeltaX_);
    ptMiss_.p4.SetPy(ptMiss_.p4.Py() - **srcUnclEnergyUpDeltaY_);
  }
  double det_cov = (*srcCOVxx_)*(*srcCOVyy_) - (*srcCOVxy_)*(*srcCOVxy_);
  double ncov_xx = (*srcCOVyy_)/det_cov;
  double ncov_xy = -(*srcCOVxy_)/det_cov;
  double ncov_yy = (*srcCOVxx_)/det_cov;
  ptMiss_.significance_corrected = ptMiss_.p4.Px() * ptMiss_.p4.Px() * ncov_xx + 2 * ptMiss_.p4.Px() * ptMiss_.p4.Py() * ncov_xy + ptMiss_.p4.Py() * ptMiss_.p4.Py() * ncov_yy;
}