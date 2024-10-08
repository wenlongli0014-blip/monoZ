#include <PileUpWeight.h>

#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TVectorD.h>

#include <HZZException.h>
#include <Logger.h>
#include <Utils.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>

PileUpWeight::PileUpWeight(
    Dataset &dataset, Options const &options, RunSampler const *runSampler)
    : cache_{dataset.Reader()}, runSampler_{runSampler},
      mu_{dataset.Reader(), "Pileup_nTrueInt"} {

  YAML::Node const config = options.GetConfig()["pileup_weight"];
  if (not config)
    throw HZZException{
        "Master configuration does not contain section \"pileup_weight\"."};

  std::filesystem::path const dataPath = Options::NodeAs<std::string>(
      config, {"data_profile"});
  LoadDataProfiles(dataPath);
  LoadSimProfile(config, dataset.Info().Name());

  // The default weight index is chosen based on the requested systematic
  // variation
  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "pileup_up")
    defaultWeightIndex_ = 1;
  else if (systLabel == "pileup_down")
    defaultWeightIndex_ = 2;
  else
    defaultWeightIndex_ = 0;
  LOG_DEBUG << "Index of default pileup weight: " << defaultWeightIndex_;
}


std::string_view PileUpWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "pileup_up";
    case 1:
      return "pileup_down";
    default:
      return "";
  }
}


void PileUpWeight::LoadDataProfiles(std::filesystem::path const &path) {
  TFile inputFile{path.c_str()};
  if (inputFile.IsZombie()) {
    HZZException exception;
    exception << "Could not open file " << path << ".";
    throw exception;
  }

  for (TObject *key : *inputFile.GetListOfKeys()) {
    auto directory = dynamic_cast<TDirectory *>(
        dynamic_cast<TKey *>(key)->ReadObj());
    if (not directory)  // The object is not a directory
      continue;

    Era era;

    std::unique_ptr<TVectorD> runRange(directory->Get<TVectorD>("run_range"));
    if (not runRange) {
      HZZException exception;
      exception << "Could not read mandatory entry \"run_range\" in directory "
          << "\"" << directory->GetName() << "\" in file " << path << ".";
      throw exception;
    }
    if (runRange->GetNoElements() != 2) {
      HZZException exception;
      exception << "Array \"run_range\" in directory \"" << directory->GetName()
          << "\" in file " << path << " must contain 2 elements but "
          << runRange->GetNoElements() << " were found.";
      throw exception;
    }
    era.minRun = std::lround((*runRange)[0]);
    era.maxRun = std::lround((*runRange)[1]);

    int profileIndex = 0;
    for (auto const &label : {"nominal", "up", "down"}) {
      auto hist = directory->Get<TH1>(label);
      if (not hist) {
        HZZException exception;
        exception << "Mandatory histogram \"" << label << "\" not found in "
            << "directory \"" << directory->GetName() << "\" in file "
            << path << ".";
        throw exception;
      }
      hist->SetDirectory(nullptr);
      era.dataProfiles[profileIndex].reset(hist);
      ++profileIndex;
    }

    eras_.emplace_back(std::move(era));
  }

  if (eras_.empty()) {
    HZZException exception;
    exception << "No data pileup profiles have been read from file "
        << path << ".";
    throw exception;
  }
  LOG_DEBUG << "Pileup profiles for " << eras_.size() << " eras loaded.";

  std::sort(
      eras_.begin(), eras_.end(),
      [](Era const &lhs, Era const &rhs){return lhs.minRun < rhs.minRun;});

  // Make sure the profiles are normalized to represent probability density
  for (auto &era : eras_)
    for (auto &profile : era.dataProfiles)
      profile->Scale(1. / profile->Integral(), "width");
}


void PileUpWeight::LoadSimProfile(
    YAML::Node const &config, std::string const &datasetName) {
  if (config["sim_profiles"]) {
    simProfile_ = utils::ReadHistogram(
        config["sim_profiles"].as<std::string>(), datasetName, false);
    if (simProfile_) {
      LOG_DEBUG << "Will use dataset-specific pileup profile in simulation.";
    } else {
      if (not config["default_sim_profile"]) {
        HZZException exception;
        exception << "File \"" << config["sim_profiles"].as<std::string>()
            << "\" does not contain pileup profile for dataset \""
            << datasetName << "\" and no default pileup profile is "
            "provided in the master configuration.";
        throw exception;
      }
      simProfile_ = utils::ReadHistogram(
          config["default_sim_profile"].as<std::string>(), "pileup");
      LOG_DEBUG << "Will use default pileup profile in simulation.";
    }
  } else {
    if (not config["default_sim_profile"])
      throw HZZException{
          "Illegal master configuration. Section \"pileup_weight\" must "
          "contain at least one of keys \"sim_profiles\" and "
          "\"default_sim_profile\"."};
    simProfile_ = utils::ReadHistogram(
        config["default_sim_profile"].as<std::string>(), "pileup");
  }

  // Normalize pileup profile to represent probability density
  simProfile_->Scale(1. / simProfile_->Integral(), "width");
}


void PileUpWeight::Update() const {
  // Note that the expected number of pileup interactions stored in NanoAOD is
  // truncated (as in std::trunc, not std::round) to an integer. Below assume
  // that profiles in data and simulation have integer binnings and make use of
  // the fact that in TH1 left boundary of a bin is inclusive. Given that small
  // integers are representable in float32 exactly, the probabilities below can
  // be computed by a simple lookup from the histograms. If this were not the
  // case, the digitization imposed by the truncation would have to be treated
  // explicitly.
  double const probSim = simProfile_->GetBinContent(
      simProfile_->FindFixBin(*mu_));
  if (probSim <= 0.) {
    LOG_WARN << "Got pileup probability in simulation of " << probSim <<
      " for true pileup of " << *mu_ << ". Set pileup weights to 1.";
    std::fill(weights_.begin(), weights_.end(), 1.);
  }

  auto const run = runSampler_->Get();
  auto const eraIter = std::lower_bound(
      eras_.begin(), eras_.end(), run,
      [](Era const &era, int const run){return era.maxRun < run;});
  if (eraIter == eras_.end() or run < eraIter->minRun) {
    HZZException exception;
    exception << "Have not found era with data pileup profiles for "
        << "run " << run << ".";
    throw exception;
  }
  auto const &era = *eraIter;

  for (int i = 0; i < int(weights_.size()); ++i) {
    double const probData = era.dataProfiles[i]->GetBinContent(
        era.dataProfiles[i]->FindFixBin(*mu_));
    weights_[i] = probData / probSim;
  }
  LOG_TRACE << "Pileup weights: mu " << *mu_ << ", run " << run << " -> "
      << weights_[0] << ", " << weights_[1] << ", " << weights_[2];
}
