#include <PileUpIdWeight.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <string>

#include <TFile.h>

#include <FileInPath.h>


/**
 * \brief Auxiliary class to simplify reading of histograms in
 * PileUpIdWeight::LoadScaleFactors
 *
 * It reads 2D histograms for a ROOT file, wrapping them in std::shared_ptr. If
 * a histogram with the same name is requested again, the already existing
 * std::shared_ptr is returned.
 */
class HistReader {
 public:
  HistReader(std::filesystem::path const &path);
  ~HistReader();
  std::shared_ptr<TH2> operator()(std::string const &name);

 private:
  TFile file_;
  std::map<std::string, std::shared_ptr<TH2>> readHistograms_;
};

HistReader::HistReader(std::filesystem::path const &path)
    : file_{path.c_str()} {}

HistReader::~HistReader() {
  file_.Close();
}

std::shared_ptr<TH2> HistReader::operator()(std::string const &name) {
  auto const res = readHistograms_.find(name);
  if (res != readHistograms_.end())
    return res->second;
  auto hist = std::shared_ptr<TH2>(file_.Get<TH2>(name.c_str()));
  hist->SetDirectory(nullptr);
  readHistograms_[name] = hist;
  return hist;
}


PileUpIdWeight::PileUpIdWeight(
    Dataset &dataset, Options const &options,
    PileUpIdFilter const *pileUpIdFilter,
    JetBuilder const *jetBuilder)
    : pileUpIdFilter_{pileUpIdFilter}, jetBuilder_{jetBuilder},
      absEtaEdges_{pileUpIdFilter_->GetAbsEtaEdges()},
      expPileUp_{dataset.Reader(), "Pileup_nTrueInt"},
      cache_{dataset.Reader()},
      systLabel_(options.GetAs<std::string>("syst")){
  for (auto const &wp : pileUpIdFilter_->GetWorkingPoints())
    contexts_.emplace_back(wp);

  std::string year = Options::NodeAs<std::string>(options.GetConfig(), {"period"});
  auto const isAPVNode = options.GetConfig()["is_APV"];
  bool const isAPV = (isAPVNode and not isAPVNode.IsNull() and isAPVNode.as<bool>());
  if (year == "2016") {
    if (isAPV) {
      year = "2016APV";
    }
  }

  bool const isUL = true; // TODO: read from config
  LoadScaleFactors(options.GetConfig(), year.c_str(), isUL);

  // auto const effModelPath = FileInPath::Resolve(
  //     Options::NodeAs<std::string>(
  //         options.GetConfig(), {"pileup_id", "efficiency"}));
  // effCalc_.emplace(effModelPath, effFeatures_.size());

  // // The year won't change, so the corresponding features can be set already
  // // now
  // effFeatures_[4] = (year == 2016) ? 1 : 0;
  // effFeatures_[5] = (year == 2017) ? 1 : 0;
  // effFeatures_[6] = (year == 2018) ? 1 : 0;
  auto const systLabel = options.GetAs<std::string>("syst");
  // if (systLabel == "puid_tag_up")
  //   defaultVariation_ = Variation::kTagUp;
  // else if (systLabel == "puid_tag_down")
  //   defaultVariation_ = Variation::kTagDown;
  // else if (systLabel == "puid_mistag_up")
  //   defaultVariation_ = Variation::kMistagUp;
  // else if (systLabel == "puid_mistag_down")
  //   defaultVariation_ = Variation::kMistagDown;
  // else
    defaultVariation_ = Variation::kNominal;
}


// std::string_view PileUpIdWeight::VariationName(int variation) const {
//   switch (variation) {
//     case 0:
//       return "puid_tag_up";
//     case 1:
//       return "puid_tag_down";
//     case 2:
//       return "puid_mistag_up";
//     case 3:
//       return "puid_mistag_down";
//     default:
//       return "";
//   }
// }


PileUpIdWeight::Context const &PileUpIdWeight::FindContext(
    Jet const &jet) const {
  int const bin = std::upper_bound(
      absEtaEdges_.begin(), absEtaEdges_.end(), std::abs(jet.p4.Eta()))
      - absEtaEdges_.begin();
  return contexts_[bin];
}


double PileUpIdWeight::GetEfficiency(
    Context const &context, Jet const &jet) const {

  std::shared_ptr<TH2> histValue;
  // std::shared_ptr<TH2> histValue, histUnc;

  histValue = context.eff;

  int const bin = histValue->FindFixBin(jet.p4.Pt(), jet.p4.Eta());
  double const effNominal = histValue->GetBinContent(bin);

  return effNominal;

}


double PileUpIdWeight::GetScaleFactor(
    Context const &context, Jet const &jet) const {
  std::shared_ptr<TH2> histValue, histValue_unc;
  // std::shared_ptr<TH2> histValue, histUnc;

  histValue = context.sf;

  int const bin = histValue->FindFixBin(jet.p4.Pt(), jet.p4.Eta());
  double sfNominal = histValue->GetBinContent(bin);

  histValue_unc = context.sf_unc;
  int const bin_unc = histValue_unc->FindFixBin(jet.p4.Pt(), jet.p4.Eta());
  double const sf_unc = histValue_unc->GetBinContent(bin_unc);

  if (systLabel_ == "puid_up")
    sfNominal += sf_unc;
  if (systLabel_ == "puid_down")
    sfNominal -= sf_unc;

  return sfNominal;
}


void PileUpIdWeight::LoadScaleFactors(YAML::Node const config, std::string year, bool isUL) {
  std::filesystem::path const path = FileInPath::Resolve(
      Options::NodeAs<std::string>(config, {"pileup_id", "scale_factors"}));
  HistReader histReader{path};
  for (auto &context : contexts_) {
    if (context.workingPoint == Jet::PileUpId::None)
      continue;

    std::string wpLabel;
    switch (context.workingPoint) {
      case Jet::PileUpId::Loose:
        wpLabel = "L";
        break;
      case Jet::PileUpId::Medium:
        wpLabel = "M";
        break;
      case Jet::PileUpId::Tight:
        wpLabel = "T";
        break;
      default:
        wpLabel = "";
    }

    std::string const nameFragment = isUL ? ("UL" + year + "_" + wpLabel) : (year + "_" + wpLabel);

    LOG_DEBUG << "Will use PileUp Jet ID eff and sf " << nameFragment << ".";
    context.sf_unc = histReader("h2_eff_sf" + nameFragment + "_Systuncty");
    context.sf = histReader("h2_eff_sf" + nameFragment);
    context.eff = histReader("h2_eff_mc" + nameFragment);

  }
}


void PileUpIdWeight::Update() const {
  std::fill(weights_.begin(), weights_.end(), 1.);

  // Jets that pass pileup ID
  for (auto const &jet : jetBuilder_->Get()) {
    if (not pileUpIdFilter_->IsTaggable(jet))
      continue;
    auto const &context = FindContext(jet);
    if (context.workingPoint == Jet::PileUpId::None)
      continue;

    // only apply to gen-matched jets
    if (jet.isPileUp)
      continue;

    double const eff = GetEfficiency(context, jet);
    double const sf = GetScaleFactor(context, jet);

    // for (int iVar = 0; iVar < 5; ++iVar) {
    for (int iVar = 0; iVar < 1; ++iVar) {
      weights_[iVar] *= std::min(sf * eff, 1.) / eff;
    }
  }

  // Jets that fail pileup ID
  for (auto const &jet : jetBuilder_->GetRejected()) {
    if (not pileUpIdFilter_->IsTaggable(jet))
      continue;
    auto const &context = FindContext(jet);
    if (context.workingPoint == Jet::PileUpId::None)
      continue;

    // only apply to gen-matched jets
    if (jet.isPileUp)
      continue;

    double const eff = GetEfficiency(context, jet);
    double const sf = GetScaleFactor(context, jet);

    // for (int iVar = 0; iVar < 5; ++iVar) {
    for (int iVar = 0; iVar < 1; ++iVar) {
      weights_[iVar] *= std::max(1. - sf * eff, 0.) / (1. - eff);
    }
  }

  LOG_TRACE << "Nominal weight in PileUpIdWeight: " << weights_[0] << ".";
}

