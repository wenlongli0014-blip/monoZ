#include <BTagWeight.h>
#include <HZZException.h>
#include <FileInPath.h>

#include <cmath>
#include <cstdlib>

#include <TFile.h>

#include "BTag/BTagCalibrationStandalone.h"

using namespace std::string_literals;


BTagWeight::BTagWeight(Dataset &dataset, Options const &options,
                       BTagger const *bTagger, JetBuilder const *jetBuilder)
    : bTagger_{bTagger}, jetBuilder_{jetBuilder},
      effTablesPath_{FileInPath::Resolve(
        Options::NodeAs<std::string>(
          options.GetConfig(), {"b_tag_weight", "efficiency"}))},
      bottomHistName_{"b"},
      charmHistName_{"c"},
      lightHistName_{"udsg"},
      scaleFactorReader_{new BTagCalibrationReader{
        // BTagEntry::OP_LOOSE, "central", {"up", "down"}}},
        BTagEntry::OP_MEDIUM, "central", {"up", "down"}}},
      cache_{dataset.Reader()} {

  std::string const scaleFactorsPath{FileInPath::Resolve(
    Options::NodeAs<std::string>(
      options.GetConfig(), {"b_tag_weight", "scale_factors"}))};

  BTagCalibration calibration{"", scaleFactorsPath};

  scaleFactorReader_->load(calibration, BTagEntry::FLAV_B, "mujets");
  scaleFactorReader_->load(calibration, BTagEntry::FLAV_C, "mujets");
  scaleFactorReader_->load(calibration, BTagEntry::FLAV_UDSG, "incl");

  auto const config = Options::NodeAs<YAML::Node>(
        options.GetConfig(), {"b_tag_weight"});

  if (auto const node = config["bottom_hist_name"]) {
    bottomHistName_ = node.as<std::string>();
  }
  if (auto const node = config["charm_hist_name"]) {
    charmHistName_ = node.as<std::string>();
  }
  if (auto const node = config["light_hist_name"]) {
    lightHistName_ = node.as<std::string>();
  }

  LoadEffTables();

  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "btag_up")
    defaultVariation_ = Variation::kTagUp;
  else if (systLabel == "btag_down")
    defaultVariation_ = Variation::kTagDown;
  else if (systLabel == "bmistag_up")
    defaultVariation_ = Variation::kMistagUp;
  else if (systLabel == "bmistag_down")
    defaultVariation_ = Variation::kMistagDown;
  else
    defaultVariation_ = Variation::kNominal;
}


BTagWeight::~BTagWeight() noexcept {}


std::string_view BTagWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "btag_up";
    case 1:
      return "btag_down";
    case 2:
      return "bmistag_up";
    case 3:
      return "bmistag_down";
    default:
      return "";
  }
}


double BTagWeight::ComputeWeight(Variation variation) const {
  double weight = 1.;

  for (auto const &jet : jetBuilder_->Get()) {
    if (not bTagger_->IsTaggable(jet))
      continue;

    double const sf = GetScaleFactor(
      jet.p4.Pt(), jet.p4.Eta(), jet.hadronFlavour, variation);

    if ((*bTagger_)(jet)) {
      weight *= sf;
    } else {
      double const eff = GetEfficiency(
          jet.p4.Pt(), jet.p4.Eta(), jet.hadronFlavour);
      weight *= (1 - sf * eff) / (1 - eff);
    }
  }

  return weight;
}


double BTagWeight::GetEfficiency(double pt, double eta, int flavour) const {

  std::string flavourLabel;
  eta = std::fabs(eta);

  switch (std::abs(flavour)) {
    case 5:
      flavourLabel = "b";
      break;

    case 4:
      flavourLabel = "c";
      break;

    default:
      flavourLabel = "udsg";
  }
  double const max_pt = effTables_.at(flavourLabel)->GetXaxis()->GetBinCenter(effTables_.at(flavourLabel)->GetNbinsX());
  int const globalBin = effTables_.at(flavourLabel)->FindFixBin(std::min(pt, max_pt), eta);
  return effTables_.at(flavourLabel)->GetBinContent(globalBin);
}


void BTagWeight::LoadEffTables() {

  TFile inputFile{effTablesPath_.c_str()};

  if (inputFile.IsZombie()) {
    HZZException exception;
    exception << "Could not open file " << effTablesPath_ << ".";
    throw exception;
  }

  for(std::string const &flavor : {"b", "c", "udsg"}) { 
    if (flavor == "b"){
      effTables_[flavor].reset(inputFile.Get<TH2>(bottomHistName_.c_str()));
    }
    else if (flavor == "c")
      effTables_[flavor].reset(inputFile.Get<TH2>(charmHistName_.c_str()));
    else if (flavor == "udsg")
      effTables_[flavor].reset(inputFile.Get<TH2>(lightHistName_.c_str()));
    
    if (not effTables_[flavor]) {
      HZZException exception;
      exception << "File " << effTablesPath_ <<
        " does not contain required histogram \"" << flavor << "\".";
      throw exception;
    }
    
    effTables_[flavor]->SetDirectory(nullptr);
  }
  inputFile.Close();
}


double BTagWeight::GetScaleFactor(double pt, double eta, int flavour,
                                  Variation variation) const {
  BTagEntry::JetFlavor translatedFlavour;
  switch (std::abs(flavour)) {
    case 5:
      translatedFlavour = BTagEntry::FLAV_B;
      break;
    case 4:
      translatedFlavour = BTagEntry::FLAV_C;
      break;
    default:
      translatedFlavour = BTagEntry::FLAV_UDSG;
  }

  std::string version{"central"};
  if (translatedFlavour != BTagEntry::FLAV_UDSG) {
    if (variation == Variation::kTagUp)
      version = "up";
    else if (variation == Variation::kTagDown)
      version = "down";
  } else {
    if (variation == Variation::kMistagUp)
      version = "up";
    else if (variation == Variation::kMistagDown)
      version = "down";
  }

  return scaleFactorReader_->eval_auto_bounds(
    version, translatedFlavour, eta, pt);
}


void BTagWeight::Update() const {
  for (int i = 0; i < int(weights_.size()); ++i)
    weights_[i] = ComputeWeight(Variation(i));
}

