#include <PileUpIdFilter.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

#include <HZZException.h>
#include <Logger.h>


PileUpIdFilter::PileUpIdFilter(Options const &options) {
  auto const configNode = options.GetConfig()["pileup_id"];
  if (not configNode) {
    enabled_ = false;
    LOG_DEBUG << "PileUpIdFilter is created but disabled because no "
        << "configuration for is has been provided.";
  } else {
    enabled_ = true;

    auto const rangeNode = Options::NodeAs<YAML::Node>(
        configNode, {"pt_range"});
    if (not rangeNode.IsSequence() or rangeNode.size() != 2)
      throw HZZException{
          "Field \"pt_range\" in section \"pileup_id\" must be a sequence of "
          "length 2."};
    minPt_ = rangeNode[0].as<double>();
    maxPt_ = rangeNode[1].as<double>();
    if (minPt_ >= maxPt_)
      throw HZZException{
          "Wrong ordering in field \"pt_range\" in section \"pileup_id\"."};

    // Field "abs_eta_edges" is optional. If not given, the same working point
    // will be used everywhere. In this case vector absEtaEdges_ is left empty.
    auto const edgesNode = configNode["abs_eta_edges"];
    if (edgesNode) {
      absEtaEdges_ = edgesNode.as<std::vector<double>>();
      if (not std::is_sorted(
            absEtaEdges_.begin(), absEtaEdges_.end(),
            [](double const &a, double const &b){return a < b;}))
        throw HZZException{
            "Sequence \"abs_eta_edges\" in section \"pileup_id\" must be "
            "sorted."};
    }

    auto const wpLabels = Options::NodeAs<std::vector<std::string>>(
        configNode, {"working_points"});
    if (wpLabels.size() != absEtaEdges_.size() + 1) {
      HZZException exception;
      exception << "In section \"pileup_id\", got " << wpLabels.size()
          << " working points and " << absEtaEdges_.size() << " |eta| edges. "
          << "The difference between the two numbers must be exactly 1.";
      throw exception;
    }
    for (auto const &label : wpLabels) {
      Jet::PileUpId wp;
      if (label == "N")
        wp = Jet::PileUpId::None;
      else if (label == "L")
        wp = Jet::PileUpId::Loose;
      else if (label == "M")
        wp = Jet::PileUpId::Medium;
      else if (label == "T")
        wp = Jet::PileUpId::Tight;
      else {
        HZZException exception;
        exception << "Illegal label \"" << label << "\" found in field "
            << "\"pileup_id\"/\"working_points\".";
        throw exception;
      }
      workingPoints_.emplace_back(wp);
    }

    LOG_DEBUG << "PileUpIdFilter is created and enabled.";
  }
}


bool PileUpIdFilter::operator()(Jet const &jet) const {
  if (not enabled_)
    return true;

  if (not (jet.p4.Pt() >= minPt_ and jet.p4.Pt() < maxPt_))
    return true;

  int const bin = std::upper_bound(
      absEtaEdges_.begin(), absEtaEdges_.end(), std::abs(jet.p4.Eta()))
      - absEtaEdges_.begin();
  return int(jet.pileUpId) >= int(workingPoints_[bin]);
}
